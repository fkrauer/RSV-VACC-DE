# Household stuff
using RSVVACCDE
using DrWatson
using Random
using DifferentialEquations, DiffEqSensitivity
using CSV, DataFrames
using Turing
using LabelledArrays
using PreallocationTools
using ForwardDiff
using Optim
using Serialization
using BenchmarkTools
using StaticArrays
using XLSX

#using ReverseDiff
#using Memoization

using StatsPlots
#using StatsBase
using Distributions
#using DistributionsAD
#using DelimitedFiles
#using BenchmarkTools
#using Optim
#using Zygote
#using StatProfilerHTML
#using ProfileView

#using LinearAlgebra
Random.seed!(1)

# Input
n_ages = 25
n_levels = 5
modelfunc = "MSEIARS5_fit" 

strategy = 0
if modelfunc == "MSEIARS5_fit"
    strategy = 0
else
    nothing
end

vₛ = true
if modelfunc == "MSEIARS5_fit"
    simyears = 4
else
    simyears = 5
end

if modelfunc == "MSEIARS5_fit"
    n_states = 10
    n_strata = 2
else
    n_states = 15
    n_strata = 6
end

# Get data
epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c  = RSVVACCDE.load_data(; n_agegroups=n_ages);

# Get params
ϵ, μ, n, pₛ, δ, qₘₑₐₙ, hₘₑₐₙ, h, agemid, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, vₙ, θ, vₛ = RSVVACCDE.load_params(epidata; strategy=strategy);

# Get inits
u0 = RSVVACCDE.load_inits(n, θ.pₕ, n_levels, n_states, n_strata, 1.0);

# get cache for preallocated version
if modelfunc == "MSEIARS5_fit" 
    cache = RSVVACCDE.make_cache_fit(n_ages, n_levels, n_strata);
else
    cache = RSVVACCDE.make_cache_vacc(n_ages, n_levels, n_strata);
end

# Define FOI forcing function
const betafunc = cosine

# Solver settings
if modelfunc == "MSEIARS5_fit" 
    tmin = minimum(times_obsindex)
    tmax = maximum(times_obsindex)
    tstart = tmin - 20.0 * 365.0
    tspan = (tstart, tmax) # fit, includes 30 years run-in for stable periodic orbit
    saveat = [(tmin-7.0):1.0:tmax;] # fit
    idxs = nothing # fit
else
    tmax = simyears*365.0
    tstart = - 20.0 * 365.0
    tspan = (tstart, tmax) # vaccsim
    saveat = [0.0:1.0:tmax;] # vaccsim
    #saveat = [0.0:365.0:tmax;] # vaccsim
    #saveat = [0.0, tmax] # vaccsim
    start_stratum = collect(range(start = 10*n_ages*n_levels + 1, step=n_ages*n_levels*n_states, length=n_strata))
    #idxs = Int.(reduce(vcat,([collect(i:1:(4*n_ages*n_levels + i -1);) for i in start_stratum])))
    idxs = nothing # fit
end

solver = Tsit5() #AutoTsit5(Rosenbrock23())

# Callback for seasonal prophylactics and/or vacc
modelfunc == "MSEIARS5_fit" ? switch = 9 : switch = 10
cb_set = load_cb_intervention(θ, 
                                iₚₐ, 
                                iₚₗ,
                                iᵥₗ,
                                iᵥₐ;
                                vₙ=vₙ,
                                vₛ = vₛ,
                                simyears=simyears,
                                modelfunc=modelfunc, 
                                switch=switch)

solvsettings = (abstol=1.0e-4,
                reltol=1.0e-4,
                save_idxs=idxs,
                saveat=saveat)
                #isoutofdomain = (u,p,t)->any(<(0),u))
                #isoutofdomain=(u,p,t)->any(<(0),u))
                #adaptive=false,
                #dt=1.0)
# 
# isoutofdomain = (u,p,t)-> any(x->x<0, u)) 

# Initiate ODE problem
if modelfunc == "MSEIARS5_fit"
    f = eval(Meta.parse(string(modelfunc, "!(; c, ϵ, μ, pₛ, iₘ, iₚₐ, iₚₗ, betafunc)")))
else
    f = eval(Meta.parse(string(modelfunc, "!(; c, ϵ, μ, pₛ, h, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, betafunc, vₙ, vₛ)")))
end
problem = ODEProblem(f, u0, tspan, (θ, δ, cache...))

# Solve
sol =  solve(problem, 
            solver; 
            callback=cb_set, 
            solvsettings...)

@assert SciMLBase.successful_retcode(sol) "failed to solve ODE for default parameters"
#=            
@benchmark solve(problem, 
            solver; 
            callback=cb_set, 
            solvsettings...);
=#

# plausibility checks

# constant number of individuals in each age group
#foo = dropdims(sum(Array(sol)[:,:,1:8,:,:], dims=(2,3,4)), dims=(2,3,4))
#moo = minimum(foo, dims=2) .- maximum(foo, dims=2)
#plot(moo)


# Small sols, fit model

# colnames
#=
statenames = ["CS", "CT", "CH", "D"]
#statenames = ["M1", "M2", "S", "E", "I", "A", "R1", "R2", "X", "CI", "CA", "CH", "D"]
agenames = ["_A"] .* string.(1:25)
levelnames = ["_L"] .* string.(0:4)
stratanames = ["_Str"] .* string.(1:6)

colnames = repeat(statenames; inner = length(agenames)*length(levelnames)) .* 
                repeat(repeat(agenames; outer = length(levelnames)) .* repeat(levelnames; inner = length(agenames)), outer = length(statenames))
colnames = repeat(colnames; outer=length(stratanames)) .* repeat(stratanames; inner = length(colnames))

df = DataFrame(sol)
rename!(df, 2:size(df,2) .=> colnames)
CSV.write("sim.csv", df)
=#

#for large sols, vacc model with all states


array = Array(sol);
indices=[Symbol("index$i") for i in 1:ndims(array)];
indexarray=[(;zip(indices,Tuple(t))...) for t in CartesianIndices(array)];
df = DataFrame(indexarray);
df.values=array[:];
CSV.write("sim.csv", df)


#=
if modelfunc == "MSEIARS5_vacc"
    statenames = ["M1", "M2", "S", "E", "I", "A", "R1", "R2", "X", "CI", "D", "CA", "CH", "T3", "SX"]
else
    statenames = ["M1", "M2", "S", "E", "I", "A", "R1", "R2", "X", "CI", "D"]
end

agenames = ["_A"] .* string.(1:25)
levelnames = ["_L"] .* string.(0:4)

if modelfunc == "MSEIARS5_vacc"
    stratanames = ["_Str"] .* string.(1:4)
else
    stratanames = ["_Str"] .* string.(1:2)
end

colnames = repeat(statenames; inner = length(agenames)*length(levelnames)) .* 
                repeat(repeat(agenames; outer = length(levelnames)) .* repeat(levelnames; inner = length(agenames)), outer = length(statenames))

colnames = repeat(colnames; outer=length(stratanames)) .* repeat(stratanames; inner = length(colnames))

rename!(df, 2:size(df,2) .=> colnames)
CSV.write("sol.csv", df)
=#

#CSV.write("sol.csv", DataFrame(sol))
#@benchmark solve(problem, solver; solvsettings...)

# Prior
prior = Core.eval(RSVVACCDE, :prior)

# Fixed theta
theta_fix = RSVVACCDE.fixedparameters(θ, prior())

# Set up as Turing model
Turing.setadbackend(:forwarddiff)
#Turing.setadbackend(:reversediff)
#Turing.setrdcache(true)

#using Logging
#Logging.disable_logging(Logging.Warn);

# Setup and condition model
model = turingmodel(prior,
                    times_obsindex,
                    ts_obsindex,
                    N_ts,
                    hosp_obsindex,
                    dim_q,
                    N_hosp,
                    qₘₑₐₙ, 
                    hₘₑₐₙ,
                    agemid,
                    N_noconv,
                    #N_neut,
                    theta_fix,
                    problem,
                    cache,
                    cb_set,
                    solver,
                    solvsettings;
                    sensealg=ForwardDiffSensitivity()) | (obs_ts = data_ts,
                                                        obs_ts_prop = data_ts_prop,
                                                        obs_hosp = data_hosp,
                                                        obs_hosp_prop = data_hosp_prop,
                                                        obs_noconv = data_noconv,
                                                        #obs_neut = data_neut,
                                                        );

#InterpolatingAdjoint(autojacvec=ZygoteVJP()),)
retval, issuccess = model();


#ProfileView.@profview retval, issuccess = model()


#=
@code_warntype model.f(
    model,
    Turing.VarInfo(model),
    Turing.SamplingContext(
        Random.GLOBAL_RNG, Turing.SampleFromPrior(), Turing.DefaultContext(),
    ),
    model.args...,
)
=#

# Prior predictive --------------------------------------------------------
model_decond = DynamicPPL.decondition(model, :obs_ts,
                                            :obs_ts_prop,
                                            :obs_hosp,
                                            :obs_hosp_prop,
                                            :obs_noconv,
                                            );

# Sample from prior
prior_sample = sample(model, Prior(), 500; progress=true)

#CSV.write("output/prior_samples.csv", DataFrame(Array(prior), :auto))

# Predict
# time series data
prior_predict_ts = quantile(predict(DynamicPPL.decondition(model, :obs_ts,), prior_sample), q=[0.025, 0.5, 0.975])
plot(prior_predict_ts[:,2]);
plot!(data_ts, color=:red);
plot!(prior_predict_ts[:,3]);
plot!(prior_predict_ts[:,4])

prior_predict_ts_prop = quantile(predict(DynamicPPL.decondition(model, :obs_ts_prop,), prior_sample), q=[0.025, 0.5, 0.975])
plot(prior_predict_ts_prop[:,2]);
plot!(data_ts_prop, color=:red);
plot!(prior_predict_ts_prop[:,3]);
plot!(prior_predict_ts_prop[:,4])

# hosp data
prior_predict_hosp = quantile(predict(DynamicPPL.decondition(model, :obs_hosp,), prior_sample), q=[0.025, 0.5, 0.975])
plot(prior_predict_hosp[:,2]);
plot!(data_hosp, color=:red);
plot!(prior_predict_hosp[:,3]);
plot!(prior_predict_hosp[:,4])

prior_predict_hosp_prop = quantile(predict(DynamicPPL.decondition(model, :obs_hosp_prop,), prior_sample), q=[0.025, 0.5, 0.975])
plot(prior_predict_hosp_prop[:,2]);
plot!(data_hosp_prop, color=:red);
plot!(prior_predict_hosp_prop[:,3]);
plot!(prior_predict_hosp_prop[:,4])

prior_predict_noconv = quantile(predict(DynamicPPL.decondition(model, :obs_noconv,), prior_sample), q=[0.025, 0.5, 0.975])
plot(prior_predict_noconv[:,2]);
plot!(data_noconv, color=:red);
plot!(prior_predict_noconv[:,3]);
plot!(prior_predict_noconv[:,4])


#=
prior_predict_hosp = Array(quantile(predict(DynamicPPL.decondition(model, :obs_hosp,), prior_sample), q=[0.025, 0.5, 0.975])[:,2:4])

# Bin according different age groups
predict_hosp_array = Array{Float64}(undef, 16*n_ages,3)
fill!(predict_hosp_array,NaN)

data_hosp_array = Array{Float64}(undef, 16*n_ages,3)
fill!(data_hosp_array,NaN)
for i in 1:size(prior_predict_hosp,1)
    for j in 1:3
        row = hosp_obs[i]
        predict_hosp_array[row,j] = prior_predict_hosp[i,j]
        data_hosp_array[row] = data_hosp[i]
    end
end

predict_hosp_array = [predict_hosp_array[i:min(i + 16 - 1, end),:] for i in 1:16:size(predict_hosp_array,1)]
data_hosp_array = [data_hosp_array[i:min(i + 16 - 1, end),:] for i in 1:16:size(data_hosp_array,1)]

plot_array = Any[] 
for i in 1:14
    p = scatter(predict_hosp_array[i][:,2],
                yerr = (predict_hosp_array[i][:,2] - predict_hosp_array[i][:,1], predict_hosp_array[i][:,3] - predict_hosp_array[i][:,2]),
                colour=:auto,
                markerstrokecolor = :auto,
                legend=false)
    scatter!(p, data_hosp_array[i], markersize=3, colour=:red, legend=false)
    push!(plot_array, p)
end
plot(plot_array..., layout=(5,3), size=(1200,1200))
=#



# Fit MLE ------------------------------------------------------------

mle = Optim.optimize(model, MAP(), Optim.Options(iterations=1_000_000, allow_f_increases=true))

mle = Optim.optimize(model, MLE(), Optim.Options(iterations=1_000_000, allow_f_increases=true))



# check fit
theta_est_mle = @LArray mle.values.array Tuple(keys(mle.f.vi.metadata))
p_mle = NamedTuple(vcat(theta_est_mle, theta_fix))
model_cond_mle = model_decond | p_mle;
retval, issuccess = model_cond_mle();

plot(retval.obs_ts, label="Modelled");
scatter!(data_ts, label="Observed", color=:red, alpha=0.4)

plot(retval.obs_ts_prop, label="Modelled");
plot!(data_ts_prop, label="Observed", color=:red, alpha=0.4)

plot(retval.obs_ts_prop, label="Modelled");
plot!(data_ts_prop, label="Observed", color=:red, alpha=0.4)

plot(retval.obs_hosp, label="Modelled");
plot!(data_hosp, label="Observed", color=:red, alpha=0.4)

plot(retval.obs_hosp_prop, label="Modelled");
plot!(data_hosp_prop, label="Observed", color=:red, alpha=0.4)

plot(retval.hh, label="Modelled");
plot!(retval.h, label="Modelled")


# FIT NUTS ------------------------------------------------------------

#chain = sample(model, MH(), 500_000, progress=true) #NUTS(100, 0.8)
chain = sample(model, NUTS(400, 0.8), 100; progress=true) 
#chain = sample(model, NUTS(100, 0.65), MCMCSerial(), 100, 2, progress=true) 
#chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

#chain = sample(model, NUTS(1000, 0.8), MCMCThreads(), 1000, 4, progress=true) 
#CSV.write("output/prior_samples.csv", DataFrame(Array(prior), :auto))
#chain = mapreduce(c -> sample(model, NUTS(1000, 0.65), 500, progress=true), chainscat, 1:4)   


# Check Fit --------------------------------------------------------------

model = deserialize("output/20-6-2023_0-15-15-828/model.jls");
chain = deserialize("output/20-6-2023_0-15-15-828/chain.jls");
chain = chain[501:end,:,:]

estparnames = Tuple(chain.name_map.parameters)
θ_fix_labels = ((k for k in keys(θ) if k ∉ estparnames)...,)
theta_fix = @LArray θ[collect(θ_fix_labels)] θ_fix_labels

n_samples = 10

predictive_model = DynamicPPL.decondition(model, 
                                        :obs_ts, 
                                        :obs_ts_prop,
                                        :obs_hosp, 
                                        :obs_hops_prop, 
                                        :obs_noconv,
                                        :obs_neut);


times_obsindex_predict = [minimum(model.args.times_obsindex):7:maximum(model.args.times_obsindex);]
ts_obsindex_predict = [1:1:length(times_obsindex_predict)*5;]
hosp_obsindex_predict = [1:1:16*25;]

# Update `predictive_model` to use `times_predict`.
predictive_model = RSVVACCDE.replace_args(predictive_model; 
                                                times_obsindex = times_obsindex_predict, 
                                                ts_obsindex = ts_obsindex_predict,
                                                hosp_obsindex = hosp_obsindex_predict);


parameters = DataFrame(MCMCChains.get_sections(chain, :parameters))

θ = NamedTuple(parameters[1,:])
cmodel = predictive_model | θ
retval, issuccess = cmodel()


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

# proportions adults with neutralising antibodies > 10 log2
predict_neut = DataFrame(retval.obs_neut, :auto)
rename!(predict_neut, 1:1 .=> ["sim"])
predict_neut[!, :season] = [1:1:4;]
predict_neut.sim = Float64.(predict_neut.sim) ./ Float64.(N_neut[:,1])

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

predict_h = DataFrame(reshape(retval.h, length(retval.h),1),:auto)
rename!(predict_h, 1:1 .=> ["sim"])
predict_h[!, :agegrp] = [1:1:25;]

predict_hh = DataFrame(reshape(retval.hh, length(retval.hh),1),:auto)
rename!(predict_hh, 1:1 .=> ["sim"])
predict_hh[!, :agegrp] = [1:1:12;]

# Vector with sol column names for exporting
statenames = ["M1", "M2", "S", "E", "I", "A", "R1", "R2", "X", "CI"]
agenames = ["_A"] .* string.(1:25)
levelnames = ["_L"] .* string.(0:4)
stratanames = ["_Str"] .* string.(1:2)
colnames = repeat(statenames; inner = length(agenames)*length(levelnames)) .* 
                repeat(repeat(agenames; outer = length(levelnames)) .* repeat(levelnames; inner = length(agenames)), outer = length(statenames))
colnames = repeat(colnames; outer=length(stratanames)) .* repeat(stratanames; inner = length(colnames))











plot(chain)
chain = chain[201:end,:,:]

plot(MCMCChains.get_sections(chain, :internals))

# Check fit
# Decondition
predictive_model = DynamicPPL.decondition(model, :obs_ts);
predictions = predict(predictive_model, chain);
predict_array = Array(predictions)

# quantiles
quantiles = mapslices(x -> quantile(x, [0.025, 0.5, 0.975]),
    predict_array; dims=1)
    quantiles = permutedims(quantiles, (2, 1))


bar(data_ts, label="Observed", color=:red, alpha=0.4)
scatter!(quantiles[:, 2];
    yerr=(quantiles[:, 2] - quantiles[:, 1], quantiles[:, 3] - quantiles[:, 2]),
    colour=:blue,
    markerstrokecolor = :blue,
    label="95% CrI")

