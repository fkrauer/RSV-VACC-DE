using ArgParse, DrWatson
using Setfield
using RSVVACCDE
using Random
using DifferentialEquations, DiffEqSensitivity
using CSV, DataFrames
using Turing, MCMCChains
using StatsBase, StatsPlots
using LabelledArrays
using PreallocationTools
using ForwardDiff
using Serialization
using Distributions
using Dates
using DynamicPPL
using AdvancedHMC.ProgressMeter


include(scriptsdir("utils.jl"))

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
end

include(scriptsdir("parse_args.jl"))
include(scriptsdir("chain_merge.jl"))

const verbose = args["verbose"]

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
chain, run_args, model = chain_merge(runs; load_model=true, ignorecommit=args["ignore-commit"])
verbose && @info "(✓) Done!"

# Load data
verbose && @info "(♻) Loading data and parameters..."
#n_ages = run_args["n_ages"]
const epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c  = RSVVACCDE.load_data(; n_agegroups=25);
const ϵ, μ, n, pₛ, δ, qₘₑₐₙ, hₘₑₐₙ, h, agemid, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, vₙ, θ = RSVVACCDE.load_params(epidata; strategy=0);
verbose && @info "(✓) Data and parameters loaded!"


# Predict 
verbose && @info "(♻) Predicting..."
# `decondition` model so that `observations` are now treated as random variables.
predictive_model = DynamicPPL.decondition(model, 
                                        :obs_ts, 
                                        :obs_ts_prop,
                                        :obs_hosp, 
                                        :obs_hops_prop, 
                                        :obs_noconv,
                                        );

# We also want to simulate for _all_ timesteps (also time points with missing data), not just those that we have data for.
times_obsindex_predict = [minimum(model.args.times_obsindex):7:maximum(model.args.times_obsindex);]
ts_obsindex_predict = [1:1:length(times_obsindex_predict)*5;]
hosp_obsindex_predict = [1:1:16*25;]

# Update `predictive_model` to use `times_predict`.
predictive_model = RSVVACCDE.replace_args(predictive_model; 
                                                times_obsindex = times_obsindex_predict, 
                                                ts_obsindex = ts_obsindex_predict,
                                                hosp_obsindex = hosp_obsindex_predict);

# Extract the parameters from `chain`.
parameters = DataFrame(MCMCChains.get_sections(chain, :parameters))

# Function to return the generated quantities to visualize the fit               
function generate_output(retval)

    ## Predict fitted quantities

    # weekly lab incidence time series (total cases)
    predict_ts = DataFrame(retval.obs_ts,:auto)
    rename!(predict_ts, 1:1 .=> ["sim"])
    predict_ts[!, :time] = times_obsindex_predict

    predict_ts_prop = DataFrame(reshape(retval.obs_ts_prop,5,1), :auto)
    rename!(predict_ts_prop, 1:1 .=> ["sim"])
    predict_ts_prop[!, :agegrp] = [1:1:5;]
    predict_ts_prop.sim = Float64.(predict_ts_prop.sim) ./ Float64(N_ts)

    # quarterly hospitalizations
    predict_hosp = DataFrame(retval.obs_hosp,:auto)
    rename!(predict_hosp, 1:1 .=> ["sim"])
    predict_hosp[!, :quarter] = [1:1:16;] 

    # proportions of hospitalizations by age group
    predict_hosp_prop = DataFrame(reshape(retval.obs_hosp_prop,25,1), :auto)
    rename!(predict_hosp_prop, 1:1 .=> ["sim"])
    predict_hosp_prop[!, :agegrp] = [1:1:25;]
    predict_hosp_prop.sim = Float64.(predict_hosp_prop.sim) ./ Float64(N_hosp)

    # proportions <1 year olds not seroconverted
    predict_noconv = DataFrame(retval.obs_noconv, :auto)
    rename!(predict_noconv, 1:1 .=> ["sim"])
    predict_noconv[!, :agegrp] = [1:1:11;]
    predict_noconv.sim = Float64.(predict_noconv.sim) ./ Float64.(N_noconv[:,1])


    ##  Other returned quantities (not fitted)

    # age-specific weekly lab incidence time series
    predict_ts_age = DataFrame(retval.inc_w_rep,:auto)
    rename!(predict_ts_age, 1:1 .=> ["sim"])
    predict_ts_age[!, :time] = repeat(times_obsindex_predict, outer=5)
    predict_ts_age[!, :agegrp] = repeat([1:1:5;], inner=length(times_obsindex_predict))

    # age-specific quarterly hospitalizations
    predict_hosp_age = DataFrame(retval.hosp_q_rep,:auto)
    rename!(predict_hosp_age, 1:1 .=> ["sim"])
    predict_hosp_age[!, :quarter] = repeat([1:1:16;], outer=25)
    predict_hosp_age[!, :agegrp] = repeat([1:1:25;], inner=16)

    # age-specific probability of hospitalization
    predict_hl = DataFrame(reshape(retval.hl, length(retval.hl),1),:auto)
    rename!(predict_hl, 1:1 .=> ["sim"])
    predict_hl[!, :agegrp] = [1:1:25;]

    predict_hh = DataFrame(reshape(retval.hh, length(retval.hh),1),:auto)
    rename!(predict_hh, 1:1 .=> ["sim"])
    predict_hh[!, :agegrp] = [1:1:12;]

    # age-specific proportion of hospitalised high level among all levels
    predict_prop_highrisk = DataFrame(retval.prop_highrisk,:auto)
    rename!(predict_prop_highrisk, 1:1 .=> ["sim"])
    predict_prop_highrisk[!, :agegrp] = [1:1:25;]

    #=
    # Full solution (all states)
    sol = DataFrame(retval.sol)
    rename!(sol, 1:1 .=> ["time"])
    # Subset: We're only interested in 1 full year of simulation
    sol = subset(sol, :time => ByRow(≥(0.0)))
    sol = subset(sol, :time => ByRow(<(365.0)))
    rename!(sol, 2:size(sol,2) .=> colnames)
    =#
        
    return (predict_ts, predict_ts_prop, predict_hosp, predict_hosp_prop, predict_noconv, 
            predict_ts_age, predict_hosp_age, predict_hl, predict_hh, predict_prop_highrisk) #, sol

end

# Simulate from posterior
output = (
        sim_ts = DataFrame(),
        sim_ts_prop = DataFrame(),
        sim_hosp = DataFrame(),
        sim_hosp_prop = DataFrame(),
        sim_noconv = DataFrame(),
        sim_ts_age = DataFrame(),
        sim_hosp_age = DataFrame(),
        sim_hl = DataFrame(),
        sim_hh = DataFrame(),
        sim_prop_highrisk = DataFrame()
        )
#sim_sol = DataFrame()

parameter_iterator = Tables.namedtupleiterator(parameters)
@showprogress "Generating..." for (i, θ) in enumerate(parameter_iterator)
    # Condition `model` on the parameters `θ``
    cmodel = predictive_model | θ
    # We wrap this in a `try-catch` because the call to `rand` for generating the
    # observations might fail if we have a particularly badly behaved `θ`.
    try
        # Sample from predictive.
        retval, issuccess = cmodel()
        if !issuccess
            @warn "iteration $i does not result in a successful evaluation of the model; skipped!"
        else
            predictions #=, sol=# = generate_output(retval)

            for (j, k) in enumerate(predictions)
                k[!, :replicate] .= i
                k[!, :chain] .= parameters.chain[i]
                append!(output[j], k)
            end
            #=
            sol[!, :replicate] .= i
            sol[!, :chain] .= parameters.chain[i]
            append!(sim_sol, sol)
            =#

        end
    catch e
        @warn "iteration $i resulted in an error when evaluating the model; skipped!"
    end
end

verbose && @info "(✓) Predictions done!"

# chain plots and diagnostics
verbose && @info "(♻) Generating plots..."
traceplot = plot(chain)
lpplot = plot(chain[:lp])

# Get diagnostics (ESS / RHAT)
ess_rhat_df = DataFrame(ess_rhat(chain))

verbose && @info "(✓) Plotting done!"

# Save output
verbose && @info "(♻) Calculating quantiles and saving output..."
save(string("output/figures/traceplot", ".png"), traceplot)
save(string("output/figures/lpplot", ".png"), lpplot)
CSV.write(string("output/ess_rhat", ".txt"), ess_rhat_df)

# Calculate median and 95% PPI from all simulated output
for (i,j) in enumerate(output)
    colsym = Symbol.(deleteat!(names(j), sort(indexin(["sim"; "replicate"; "chain"], names(j)))))
    out = combine(groupby(j, colsym), 
    [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95", cols*"_median", cols*"_up95"]
        for cols in ["sim"]])
    df = keys(output)[i]
    CSV.write(string("output/fit_simulations/", "quantiles_", df, ".csv"), out)
end
#CSV.write(string("output/", "sol", ".csv"), sim_sol) # this file is excessively large, do not save unless absolutely necessary
CSV.write(string("output/quantiles_posterior", ".csv"), quantile(chain))

# Draw a sample of 500 from the posterior samples and save as csv
#chaindf = DataFrame(chain)
#nsamples_chain = min(nrow(chaindf), 500)
#chainsample = chaindf[sample(1:nrow(chaindf), nsamples_chain, replace=false), :]
#CSV.write(string("output/chain_sample_500", ".csv"), chainsample)

verbose && @info "(✓) Finito!"
