using ArgParse, Dates, DrWatson

_available_metric_types = ["DiagEuclideanMetric", "DenseEuclideanMetric"]

function defaultname()
    d = Dates.now()
    return Dates.format(d, dateformat"d-m-Y_H-M-S-s")
end

s = ArgParseSettings()
@add_arg_table! s begin
    "seed"
    help = "random seed to use"
    arg_type = Int
    required = true
    "--betafunc"
    default = "cosine"
    # Solver parameters
    "--abstol"
    help = "absolute tolerance to use for `solve`"
    arg_type = Float64
    default = 1.0e-4
    "--reltol"
    help = "relative tolerance to use for `solve`"
    arg_type = Float64
    default = 1.0e-4
    "--solver"
    help = "solver to use"
    default = "Tsit5()" #"AutoTsit5(Rosenbrock23())" 
    "--chunksize"
    help = "chunk size for the automatic differentiation (default of 0 means Forwarddiff chooses the chunk size automatically)"
    arg_type = Int
    default = 0
    # Sampler parameters
    "--nadapts"
    arg_type = Int
    default = 500
    "--nsamples"
    arg_type = Int
    default = 1000
    "--max-depth"
    arg_type = Int
    default = 10
    "--target-acceptance"
    arg_type = Float64
    default = 0.8
    "--metricT"
    help = "metric from AdvancedHMC.jl to use in HMC; must be one of " * join(_available_metric_types, ", ", " or ")
    default = "DenseEuclideanMetric" #"DiagEuclideanMetric"
    range_tester = ∈(_available_metric_types)
    "--stepsize"
    help = "stepsize used for integrator in HMC; if ≤ 0.0, will use `AdvancedHMC.find_good_stepsize` to initialize the stepsize"
    default = -1.0
    arg_type = Float64
    "--integrator"
    help = "integrator from AdvancedHMC.jl to use in HMC"
    default = "AdvancedHMC.Leapfrog(stepsize)" # "AdvancedHMC.TemperedLeapfrog(stepsize, temperrate)"
    "--momentum-refreshment"
    help = "momentum refreshment from AdvancedHMC.jl to use in HMC"
    default = "AdvancedHMC.FullMomentumRefreshment()"
    "--termination-criterion"
    help = "termination criterion from AdvancedHMC.jl to use in HMC"
    default = "AdvancedHMC.GeneralisedNoUTurn(; max_depth)"
    "--stepsize-fixed"
    help = "if specified, the stepsize will not be adapted"
    action = :store_true
    "--metric-fixed"
    help = "if specified, the metric will not be adapted"
    action = :store_true
    # Model parameters
    "--initialize"
    help = "if specified, the parameters will be initialized by drawing from the prior instead of uniform"
    action = :store_true
    # Storing output
    "--outdir"
    default = "output"
    "--name"
    default = "$(defaultname())"
    # General
    "--verbose"
    help = "if specified, additional information will be logged"
    action = :store_true
    "--no-progress"
    help = "if specified, no progress meter will be shown"
    action = :store_true
end

# This will parse args if necessary, ensuring that
# a `args` variable is available.
include(scriptsdir("parse_args.jl"))

verbose = args["verbose"]
verbose && @info "args" args

seed = args["seed"]
_outdir = args["outdir"]

name = args["name"]
outdir(args...) = projectdir(_outdir, name, args...)

# Save all the inputs so that it's fully reproducible.
@tagsave outdir("args.jld2") args

# Packages
using RSVVACCDE
using Random
using DifferentialEquations, DiffEqSensitivity
using CSV, DataFrames
using Turing
using LabelledArrays
using PreallocationTools
using ForwardDiff
using Optim
using Serialization
using Distributions
using XLSX

using Turing.DynamicPPL
using AdvancedHMC
using AdvancedHMC: Adaptation

Random.seed!(seed);

const nadapts = args["nadapts"]
const nsamples = args["nsamples"]

# Prep --------------------------------------------------------------------

# Model structure
const n_ages = 25
const n_levels = 5
const n_states = 11
const n_strata = 2
const modelfunc = "MSEIARS5_fit"

# Get data
const epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c  = RSVVACCDE.load_data(; n_agegroups=n_ages);

# Get params
const ϵ, μ, n, pₛ, δ, qₘₑₐₙ, hₘₑₐₙ, h, agemid, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, vₙ, θ = RSVVACCDE.load_params(epidata; strategy=0);

# Get inits
u0 = RSVVACCDE.load_inits(n, θ.pₕ, n_levels, n_states, n_strata, 100.0);

# Set up ODE Model function
const betafunc = eval(Meta.parse(args["betafunc"]))
const f = eval(Meta.parse(string(modelfunc, "!(; c, ϵ, μ, pₛ, iₘ, iₚₐ, iₚₗ, betafunc)")))
const cache = RSVVACCDE.make_cache_fit(n_ages, n_levels, n_strata)

# Solver settings
const solver = eval(Meta.parse(args["solver"]))

tmin = minimum(times_obsindex)
tmax = maximum(times_obsindex)
tstart = tmin - 20.0 * 365.0
const tspan = (tstart, tmax) # includes 20 years run-in for stable periodic orbit
const saveat = [(tmin-7.0):1.0:tmax;] 
const solvsettings = (abstol=args["abstol"],
                        reltol=args["reltol"],
                        #isoutofdomain=(u,p,t)->any(<(0),u),
                        saveat=saveat
                        )

# Callback for seasonal prophylactics
const cb_set = load_cb_intervention(θ, 
                                    iₚₐ, 
                                    iₚₗ,
                                    iᵥₗ,
                                    iᵥₐ;
                                    vₙ=vₙ,
                                    simyears=4,
                                    modelfunc=modelfunc,
                                    switch=9)


# initialise ODE problem
const problem = ODEProblem(f, u0, tspan, (θ, δ, cache...))

# Prior
const prior = Core.eval(RSVVACCDE, :prior)

# Get Fixed theta
const theta_fix = RSVVACCDE.fixedparameters(θ, prior())

# Test ODE model
sol_test = solve(problem, solver; solvsettings...)
@assert sol_test.retcode == :Success "failed to solve ODE for default parameters" #SciMLBase.successful_retcode(sol_test) "failed to solve ODE for default parameters"

# Set up Turing model
@info "(♻) Instantiating Turingmodel..."
Turing.setadbackend(:forwarddiff)
if args["chunksize"] > 0
    Turing.setchunksize(args["chunksize"])
end

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
                                                            );

# Execute once to ensure that it's working correctly.
retval, issuccess = model();
#@assert model().obs_ts == data_ts "conditioned observations are not same as data"
@info "(✓) Turingmodel successfully instantiated!"

# Determine initial parameters.

# Author of this section: Tor Erlend Fjelde

# This ends up being quite "ugly" because we need to ensure that we're using
# the same initialization scheme as HMC for the variables that we're not going to
# explicitly set.
# 
# The way it's done is:
# 1. Initialize `VarInfo` using `SampleFromPrior()`.
# 2. `link` to move to unconstrained space.
# 3. Sample parameters in unconstrained space using `SampleFromUniform`.
# 4. `invlink!` to get back to constrained space.
prior_sampler = DynamicPPL.SampleFromPrior()
var_info = DynamicPPL.VarInfo(model, prior_sampler);

function initialize_varinfo!!(
    model,
    var_info,
    sampler=SampleFromPrior();
    initialize=args["initialize"]
)
    # `SampleFromUniform` for some reason implicitly does `invlink!`
    # even though the original `var_info` is linked. Why? Who knows.
    # But because of this, we need to keep track of if `var_info` `was_linked`
    # so that if it indeed already was linked, we need to `link!` again at the end.
    was_linked = DynamicPPL.islinked(var_info, sampler)

    # Sample initial values using `SampleUniform`.
    model(var_info, DynamicPPL.SampleFromUniform())

    if initialize
        @info "(♻) Sampling initial parameter states from prior..."
        init_params = var_info[DynamicPPL.SampleFromPrior()]
        serialize(outdir("theta_inits.jls"), init_params)
        DynamicPPL.setval!(var_info, NamedTuple{keys(var_info.metadata)}(init_params))
        @info "(✓) Done!"
    end

    # Transform back to unconstrained space if original input was unconstrained.
    was_linked && DynamicPPL.link!(var_info, sampler)

    return var_info
end


##################################
# AdvancedHMC.jl sample iterator #
##################################

# Makes it so we can use the samples from AHMC as we would a chain obtained from Turing.jl.
function Turing.Inference.HMCTransition(
    ::AdvancedHMC.HMCSampler,
    transition::AdvancedHMC.Transition,
    vi::DynamicPPL.VarInfo,
    )
    sampler = SampleFromPrior()

    # If it was not linked, we first need to convert to real space.
    !(DynamicPPL.islinked(vi, sampler)) && DynamicPPL.link!(vi, sampler)

    # Set the values in `vi`.
    var_info[sampler] = transition.z.θ

    # If it was linked, we `invlink!` to ensure that `vi` is left unchanged.
    DynamicPPL.invlink!(vi, sampler)

    # Extract `θ` and logjoint.
    θ = var_info[sampler]

    return Turing.Inference.HMCTransition(θ, transition.z.ℓπ.value, transition.stat)
end

# Need to convert the results into a `MCMCChains.Chains` object.
# We do this by overloading `AbstractMCMC.bundle_samples`.
function AbstractMCMC.bundle_samples(
    ts::Vector{<:Turing.Inference.HMCTransition{<:AbstractVector}},
    model::DynamicPPL.Model,
    chain_type::Type{MCMCChains.Chains};
    save_state=false,
    kwargs...
    )
    # Execute `model` once to get trace structure.
    var_info = DynamicPPL.VarInfo(model)

    # Convert transitions to array format.
    # Also retrieve the variable names.
    nms, _ = Turing.Inference._params_to_array([var_info])

    # Get the values of the extra parameters in each transition.
    extra_params, extra_values = Turing.Inference.get_transition_extras(ts)

    # We make our own `vals`
    vals = map(ts) do t
        Matrix(t.θ')
    end
    vals = vcat(vals...)

    # Extract names & construct param array.
    nms = [nms; extra_params]
    parray = hcat(vals, extra_values)

    info = NamedTuple()

    # Conretize the array before giving it to MCMCChains.
    parray = convert(Array{Float64}, parray)

    # Chain construction.
    return MCMCChains.Chains(
        parray,
        nms,
        (internals=extra_params,);
        info=info
    )
end

# TODO: Make PR to Turing.jl to address this.
function gen_∂logπ∂θ(vi, spl::DynamicPPL.AbstractSampler, model)
    #function Turing.Inference.gen_∂logπ∂θ(vi, spl::DynamicPPL.AbstractSampler, model)
    _vi = deepcopy(vi)
    function ∂logπ∂θ(x)
        return Turing.Core.gradient_logp(x, _vi, model, spl)
    end
    return ∂logπ∂θ
end

# TODO: Make PR to AdvancedHMC.jl to address this.
Adaptation.reset!(adaptor::Adaptation.StepSizeAdaptor) = nothing
function Adaptation.adapt!(
    ::Adaptation.StepSizeAdaptor,
    θ::AbstractVecOrMat{<:AbstractFloat},
    α::AdvancedHMC.AbstractScalarOrVec{<:AbstractFloat}
    )
    return nothing
end

# Set up the `DynamicPPL.Sampler` and `AdvancedHMC.DifferentiableDensity`.
!DynamicPPL.islinked(var_info, prior_sampler) && DynamicPPL.link!(var_info, prior_sampler)

# Sample the initial parameters.
var_info = initialize_varinfo!!(model, var_info)

# Extract the sampled parameter.
θ_init = var_info[prior_sampler]

# NOTE: We need to make sure we're in the real space when calling `gen_logπ` and others
# since these methods will return functions which have used the passed `var_info` as a
# template.
if !DynamicPPL.islinked(var_info, prior_sampler)
    error("Something went wrong; we're not in unconstrained space!")
end

# Create the logjoint and ∇logjoint computations.
logπ = Turing.Inference.gen_logπ(var_info, prior_sampler, model)
#∂logπ∂θ = Turing.Inference.gen_∂logπ∂θ(var_info, prior_sampler, model)
∂logπ∂θ = gen_∂logπ∂θ(var_info, prior_sampler, model)

while !isfinite(logπ(θ_init)) || !(all(isfinite.(∂logπ∂θ(θ_init))))
    @warn "Initial parameters lead to numerical issues; sampling some new ones!"
    global var_info, θ_init
    var_info = initialize_varinfo!!(model, var_info)
    θ_init = var_info[prior_sampler]
end

serialize(outdir("theta_inits_transformed.jls"), var_info)
# Inits in the constraint original space can be obtained with:
# DynamicPPL.invlink!(var_info, SampleFromPrior())
# DynamicPPL.getval(var_info, keys(var_info))

if verbose
    DynamicPPL.invlink!(var_info, prior_sampler)
    @info "Initial parameters (constrained)" var_info[prior_sampler]
    DynamicPPL.link!(var_info, prior_sampler)
end

# Create a model compatible with `AdvancedHMC.jl`
model_ahmc = DifferentiableDensityModel(logπ, ∂logπ∂θ)

# Construct the initial metric.
metricT = eval(Meta.parse(args["metricT"]))
metric = metricT(length(θ_init))

stepsize = if args["stepsize"] > 0 || args["stepsize-fixed"]
    args["stepsize"]
else
    verbose && @info "`stepsize` ≤ 0 → finding good initial stepsize..."
    find_good_stepsize(
        Hamiltonian(metric, model_ahmc.ℓπ, model_ahmc.∂ℓπ∂θ),
        θ_init
    )
end
verbose && @info "Using stepsize $(stepsize)"

# Evaluate sampler arguments.
const max_depth = args["max-depth"]
const target_acceptance = args["target-acceptance"]

fixed_stepsize = args["stepsize-fixed"]
fixed_metric = args["metric-fixed"]

integrator = eval(Meta.parse(args["integrator"]))
momentum_refreshment = eval(Meta.parse(args["momentum-refreshment"]))
termination_criterion = eval(Meta.parse(args["termination-criterion"]))

# Construct the initial HMC kernel.
trajectory = Trajectory{MultinomialTS}(integrator, termination_criterion)
kernel = HMCKernel(trajectory)

verbose && @info "Resulting kernel: $(kernel)"

# Construct the adaptor.
stepsize_adaptor = if fixed_stepsize
    Adaptation.FixedStepSize(AdvancedHMC.step_size(integrator))
else
    StepSizeAdaptor(target_acceptance, kernel.τ.integrator)
end
metric_adaptor = if fixed_metric
    Adaptation.UnitMassMatrix()
else
    MassMatrixAdaptor(metric)
end

adaptor = StanHMCAdaptor(
    metric_adaptor,
    stepsize_adaptor
)

# Finally construct the full sampler.
sampler_ahmc = AdvancedHMC.HMCSampler(kernel, metric, adaptor)


# Construct tempered sampler.
# No tempering, i.e always return 1.
β = i -> 1

sampler_tempered = RSVVACCDE.SerialTempering(sampler_ahmc, β)

# Iterator which we will step through.
iterator = AbstractMCMC.steps(
    model_ahmc, sampler_tempered;
    init_params=θ_init,
    nadapts=nadapts
);

# Initial step.
transition, state = iterate(iterator)

# Convert it into something we can interact with in a nicer way.
transition_turing = Turing.Inference.HMCTransition(sampler_ahmc, transition, var_info);

# Create empty sample container.
samples = AbstractMCMC.samples(
    transition_turing,
    model_ahmc,
    sampler_ahmc
);

# Save sample.
samples = AbstractMCMC.save!!(samples, transition_turing, 1, model_ahmc, sampler_ahmc);

# Set up callbacks. We want to use the `TensorBoardLogger` one + `HMCProgressCallback`.
callback_progress = if args["no-progress"]
    (args...; kwargs...) -> nothing
else
    AdvancedHMC.HMCProgressCallback(nadapts + nsamples)
end

# Execute callbacks.
callback_progress(iterator.rng, model_ahmc, sampler_ahmc, transition, state.state, length(samples); iterator.kwargs...);

verbose && @info "Adapting..."
while length(samples) < nadapts
    global transition, transition_turing, state, samples
    transition, state = iterate(iterator, state)

    # Convert it into something we can interact with in a nicer way.
    transition_turing = Turing.Inference.HMCTransition(sampler_ahmc, transition, var_info)

    # Save sample.
    samples = AbstractMCMC.save!!(samples, transition_turing, length(samples), model_ahmc, sampler_ahmc)

    # Execute callbacks.
    callback_progress(iterator.rng, model_ahmc, sampler_ahmc, transition, state.state, length(samples); iterator.kwargs...)
end

# Save the adapted state in case we want to re-use in some other run.
verbose && @info "Saving adapted state..."
serialize(outdir("state_adapt.jls"), state)

# Finally sample!
verbose && @info "Sampling..."
while length(samples) < nadapts + nsamples
    global transition, transition_turing, state, samples
    transition, state = iterate(iterator, state)

    # Convert it into something we can interact with in a nicer way.
    transition_turing = Turing.Inference.HMCTransition(sampler_ahmc, transition, var_info)

    # Save sample.
    samples = AbstractMCMC.save!!(samples, transition_turing, length(samples), model_ahmc, sampler_ahmc)

    # Execute callbacks.
    callback_progress(iterator.rng, model_ahmc, sampler_ahmc, transition, state.state, length(samples); iterator.kwargs...)
end

# Convert the `Vector` of samples into `MCMCChains.Chains`.
chain = AbstractMCMC.bundle_samples(samples, model, MCMCChains.Chains)
show(stdout, MIME"text/plain"(), chain)

# Store the results.
serialize(outdir("chain.jls"), chain)
serialize(outdir("model.jls"), model)

# Save some summarizing information allowing us to quickly
# filter later when we do `collect_results`.
ess_min = minimum(summarize(chain[nadapts+1:end]).nt.ess)
DrWatson.safesave(outdir("result.jld2"), Dict("ess_min" => ess_min))
