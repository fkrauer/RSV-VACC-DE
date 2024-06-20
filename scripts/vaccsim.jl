using ArgParse, DrWatson, ProgressBars 
using Random
using RSVVACCDE
using CSV, DataFrames, StatsPlots, DifferentialEquations, Turing, DiffEqSensitivity
using Distributions, StatsBase, LabelledArrays, MCMCChains, Serialization
using Random
using Dates
using DynamicPPL
using StaticArrays


include(scriptsdir("utils.jl"))

_available_sol_types = ["full", "quantiles"]

s = ArgParseSettings()
@add_arg_table! s begin
    "runs"
    help = "runs to process"
    action = :store_arg
    nargs = '*'
    "--prefix"
    help = "Prefix that name of run should match, e.g. '2021-08-19T14' for a runs started on 19th of August 2021 between 2pm and 3pm. Either this or runs should be provided."
    "--list-runs"
    help = "If specified, available runs will be listed."
    action = :store_true
    "--verbose"
    help = "If specified, more logging information will be provided."
    action = :store_true
    "--ignore-commit"
    help = "If specified, checkout of the correct commit won't be forced."
    action = :store_true
    "--simyears"
    help = "Duration of simulated vaccination in years"
    arg_type = Int64
    default = 5
    "--n_samples"
    help = "Number of samples to draw from posterior (for trajectory uncertainty)"
    arg_type = Int64
    default = 100
    "--strategies"
    help = "Alternative strategies to simulate"
    arg_type = String
    default = "[1, 2, 3, 4, 5, 6]"
    "--comparator"
    help = "Comparator strategy"
    arg_type = String
    default = "0"
    "--solution"
    help = "Returns either the full solution or the quantiles: " * join(_available_sol_types, ", ", " or ")
    default = "full"
    range_tester = ∈(_available_sol_types)
    "--sensitivity"
    help = "Argument to run the scenarios for the sensitivity analyses instead of the main scenarios"
    action = :store_true

end

include(scriptsdir("parse_args.jl"))
include(scriptsdir("chain_merge.jl"))

const verbose = args["verbose"]

solution = args["solution"]
args["sensitivity"] ? sensitivity = true : sensitivity = false
strategies = eval(Meta.parse(args["strategies"]))
comparator = eval(Meta.parse(args["comparator"]))

Random.seed!(42);

# Get the available runs. Useful for either `--list-runs` or if `--prefix` is given.
available_runs = filter(
    isdir ∘ Base.Fix1(projectdir, "output"),
    readdir(projectdir("output"))
)

# List runs and exit, if specified.
if args["list-runs"]
    println("The following runs are available:")
    println()
    for run in available_runs
        println("  $(run)")
    end
    exit(0)
end

# Get the absolute paths of the runs.
runs = if length(args["runs"]) > 0
    map(abspath, args["runs"])
elseif !isnothing(args["prefix"])
    prefix = args["prefix"]
    map(Base.Fix1(projectdir, "output"), filter(startswith(prefix), available_runs))
else
    error("either runs or prefix has to be provided")
end
verbose && @info "Using runs $(runs)"
@assert length(runs) > 0 "no runs specified"

# Load model, chain and model args
verbose && @info "(♻) Loading model and chain(s)..."
chain, run_args = chain_merge(runs; ignorecommit=args["ignore-commit"])
verbose && @info "(✓) Done!"

# Get data and parameters
verbose && @info "(♻) Setting up the model..."

const n_samples = min(size(chain)[1], args["n_samples"])
@assert n_samples <= 500.0 "Do not draw more than 500 samples, otherwise the resulting vaccination simulation dataframe becomes too large"
const modelfunc = "MSEIARS5_vacc"
const n_ages = 25
const n_levels = 5
const n_states = 15
const n_strata = 6
const epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c  = RSVVACCDE.load_data(; n_agegroups=n_ages);
ϵ, μ, n, pₛ, δ, qₘₑₐₙ, hₘₑₐₙ, h, agemid, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, vₙ, θ, vₛ = RSVVACCDE.load_params(epidata; strategy=0, sensitivity=sensitivity);
const u0 = RSVVACCDE.load_inits(n, θ.pₕ, n_levels, n_states, n_strata, 100.0);
const cache = RSVVACCDE.make_cache_vacc(n_ages, n_levels, n_strata);
const betafunc = cosine


# Define simulation and solver settings
simyears = args["simyears"] 
tmin = 0.0
tstart = tmin - 20.0 * 365.0
tmax = simyears*365.0
tspan = (tstart, tmax) 
saveat = [tmin:365.0:tmax;]
start_stratum = collect(range(start = 10*n_ages*n_levels + 1, step=n_ages*n_levels*n_states, length=n_strata))
idxs = Int.(reduce(vcat,([collect(i:1:(5*n_ages*n_levels + i -1);) for i in start_stratum]))) # we only save the incidence states, and the counter states for n doses and n immunised
solver = Tsit5() #AutoTsit5(Rosenbrock23()) 
solvsettings = (abstol=1.0e-4,
                reltol=1.0e-4,
                save_idxs=idxs,
                saveat=saveat)

# make vector of fixed parameters            
#theta_fix = RSVVACCDE.fixedparameters(θ, model.args.prior())
estparnames = Tuple(chain.name_map.parameters)
theta_fix_labels = ((k for k in keys(θ) if k ∉ estparnames)...,)
theta_fix = @LArray θ[collect(theta_fix_labels)] theta_fix_labels

# Set up the model
f = MSEIARS5_vacc!(; c, ϵ, μ, pₛ, h, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, betafunc, vₙ, vₛ)
problem = ODEProblem(f, u0, tspan, (θ, δ, cache...))
cb_set = load_cb_intervention(θ, 
                                iₚₐ, 
                                iₚₗ,
                                iᵥₗ,
                                iᵥₐ;
                                vₙ = vₙ,
                                vₛ = vₛ,
                                simyears=simyears,
                                modelfunc=modelfunc,
                                switch=10)
verbose && @info "(✓) Done!"


# Draw from posteriors (save only once for main analysis) --------------------------------------------------
verbose && @info "(♻) Sampling from posterior..."
posteriors = Array(chain, [:parameters])
posterior_sample = posteriors[sort(sample(1:size(posteriors, 1), n_samples, replace=false)), :]
posterior_sample_df = DataFrame(posterior_sample, :auto)
rename!(posterior_sample_df, reduce(vcat, estparnames))
if !sensitivity
    CSV.write("output/posterior_sample_100.csv", posterior_sample_df)
end
verbose && @info "(✓) Done!"


# Simulate  --------------------------------------------------------

# simulate the base scenario

verbose && @info "(♻) Simulating base (comparator) strategy ..."
inc_base, doses_base, immunisations_base = ODEwrap_vaccsim(modelfunc,
                                                            betafunc,
                                                            agemid,
                                                            hₘₑₐₙ,
                                                            c,
                                                            ϵ, 
                                                            μ, 
                                                            pₛ, 
                                                            iₘ,
                                                            u0,
                                                            tspan,
                                                            cache,
                                                            estparnames, 
                                                            posterior_sample, 
                                                            epidata, 
                                                            solvsettings, 
                                                            solver, 
                                                            n_ages, 
                                                            n_levels, 
                                                            n_strata;
                                                            strategy=comparator, 
                                                            simyears=simyears,
                                                            n_samples=n_samples,
                                                            sensitivity = sensitivity)
                                        
select!(inc_base, Not([:strategy, :seasonal]))
select!(doses_base, Not([:strategy, :seasonal]))
select!(immunisations_base, Not([:strategy, :seasonal]))
verbose && @info "(✓) Done!"


# simulate the different mAb/vaccination strategies
verbose && @info "(♻) Simulating alternative mAb and vaccination strategies ..."
inc_final = DataFrame()
doses_final = DataFrame()
immunisations_final = DataFrame()


for i in strategies
    
    verbose && @info "Simulating strategy $(i)..."

    inc, doses, immunisations = ODEwrap_vaccsim(modelfunc,
                                                betafunc,
                                                agemid,
                                                hₘₑₐₙ,
                                                c,
                                                ϵ, 
                                                μ, 
                                                pₛ, 
                                                iₘ,
                                                u0,
                                                tspan,
                                                cache,
                                                estparnames, 
                                                posterior_sample, 
                                                epidata, 
                                                solvsettings, 
                                                solver, 
                                                n_ages, 
                                                n_levels, 
                                                n_strata;
                                                strategy=i, 
                                                simyears=simyears,
                                                n_samples=n_samples,
                                                sensitivity=sensitivity)

                                
    merged_inc = outerjoin(inc, inc_base, on = [:age, :incidence, :year, :replicate], renamecols = "" => "_base") 
    merged_inc[!, :n_avert] = merged_inc.n_base .- merged_inc.n
    merged_inc[!, :prop_avert] = merged_inc.n_avert ./ merged_inc.n_base
    append!(inc_final, merged_inc) 

    merged_doses = outerjoin(doses, doses_base, on = [:age, :intervention, :year, :replicate], renamecols = "" => "_base")
    append!(doses_final, merged_doses) 

    merged_immunisations = outerjoin(immunisations, immunisations_base, on=[:age, :intervention, :year, :replicate], renamecols="" => "_base")
    append!(immunisations_final, merged_immunisations)

    
    verbose && @info "(✓) Strategy $(i) done!"

end

verbose && @info "(✓) All simulations done!"


# Calculate median and 95% CI from all sampled trajectories
if solution == "quantiles" 

    verbose && @info "(♻) Calculating quantiles of simulated trajectories ..."

    quantiles_inc = combine(groupby(inc_final, [:age, :incidence, :year, :strategy, :vaccine, :seasonal]), 
            [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95CI", cols*"_median", cols*"_up95CI"]
                for cols in ["n", "n_base", "n_avert", "prop_avert"]])

    quantiles_doses = combine(groupby(doses_final, [:age, :intervention, :year, :strategy, :vaccine, :seasonal]), 
            [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95CI", cols*"_median", cols*"_up95CI"]
                for cols in ["n", "n_base"]])

    quantiles_immunisations = combine(groupby(immunisations_final, [:age, :intervention, :year, :strategy, :vaccine, :seasonal]),
        [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols * "_low95CI", cols * "_median", cols * "_up95CI"]
        for cols in ["n", "n_base"]])

    verbose && @info "(✓) Done!"

end

verbose && @info "(♻) Saving results ..."

filename = string(comparator) * "_vs_" * string(strategies)

if sensitivity
    if solution == "full"
        CSV.write("output/vacc_simulations/vaccsim_inc_sensitivity_" * filename * ".csv", inc_final)
        CSV.write("output/vacc_simulations/vaccsim_doses_sensitivity_" * filename * ".csv", doses_final)
        CSV.write("output/vacc_simulations/vaccsim_immunisations_sensitivity_" * filename * ".csv", immunisations_final)
    else
        CSV.write("output/vacc_simulations/vaccsim_quantiles_inc_sensitivity_" * filename * ".csv", quantiles_inc)
        CSV.write("output/vacc_simulations/vaccsim_quantiles_doses_sensitivity_" * filename * ".csv", quantiles_doses)
        CSV.write("output/vacc_simulations/vaccsim_quantiles_immunisations_sensitivity_" * filename * ".csv", quantiles_immunisations)
    end
else
    if solution == "full"
        CSV.write("output/vacc_simulations/vaccsim_inc_" * filename * ".csv", inc_final)
        CSV.write("output/vacc_simulations/vaccsim_doses_" * filename * ".csv", doses_final)
        CSV.write("output/vacc_simulations/vaccsim_immunisations_" * filename * ".csv", immunisations_final)
    else
        CSV.write("output/vacc_simulations/vaccsim_quantiles_inc_" * filename * ".csv", quantiles_inc)
        CSV.write("output/vacc_simulations/vaccsim_quantiles_doses_" * filename * ".csv", quantiles_doses)
        CSV.write("output/vacc_simulations/vaccsim_quantiles_immunisations_" * filename * ".csv", quantiles_immunisations)
    end
end

verbose && @info "(✓) Done!"


