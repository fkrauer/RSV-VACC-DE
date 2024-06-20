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

include(scriptsdir("utils.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "runs"
    help = "runs to visualize"
    action = :store_arg
    nargs = '*'
    "--prefix"
    help = "Prefix that name of run should match, e.g. '2021-08-19T14' for a runs started on 19th of August 2021 between 2pm and 3pm. Either this or runs should be provided."
    # Plot/figure related options.
    "--filename"
    help = "File name of resulting visualization."
    default = "summary.pdf"
    "--theme"
    help = "Plotting theme to use."
    default = :default
    arg_type = Symbol
    # Alternatives.
    "--list-runs"
    help = "If specified, available runs will be listed."
    action = :store_true
    # General.
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
verbose && @info "(♻) Loading data..."
#n_ages = run_args["n_ages"]
const epidata, data_ts, data_ts_age, data_ts_prop, N_ts, times_obsindex, ts_obsindex, data_hosp, data_hosp_age, data_hosp_prop, N_hosp, hosp_obsindex, dim_q, data_noconv, N_noconv, c  = RSVVACCDE.load_data(; n_agegroups=25);
verbose && @info "(✓) Data loaded!"

const filename = args["filename"]

# Quantiles plotting function

# Time series data (line plot with ribbons)
plot_quantiles_ts(args...; kwargs...) = plot_quantiles_ts!(plot(), args...; kwargs...)
function plot_quantiles_ts!(p::Plots.Plot, xs, q=0.95; kwargs...)
    Δ = (1 - q) / 2
    quantiles = mapreduce(hcat, xs) do x
        quantile(x, [Δ, 0.5, 1 - Δ])
    end

    plot!(
        p,
        quantiles[2, :];
        ribbon=(quantiles[2, :] - quantiles[1, :], quantiles[3, :] - quantiles[2, :]),
        label="$(q * 100)% credible intervals",
        kwargs...
    )
    return p
end

# Rest of the data (point estimates with error bars)
plot_quantiles(args...; kwargs...) = plot_quantiles!(scatter(), args...; kwargs...)
function plot_quantiles!(p::Plots.Plot, xs, q=0.95; kwargs...)
    Δ = (1 - q) / 2
    quantiles = mapreduce(hcat, xs) do x
        quantile(x, [Δ, 0.5, 1 - Δ])
    end

    scatter!(
        p,
        quantiles[2, :];
        yerr=(quantiles[2, :] - quantiles[1, :], quantiles[3, :] - quantiles[2, :]),
        #label="$(q * 100)% credible intervals",
        markersize=0.1, # :auto,
        kwargs...
    )
    return p
end


# Display the summary of the chain.
show(stdout, MIME"text/plain"(), chain)


# chain plots.
p1 = plot(chain)

internals = filter(
    !=(:max_hamiltonian_energy_error),
    chain.name_map.internals
)
p_internals = plot(MCMCChains.get_sections(chain, :internals)[internals])

# Posterior predictive
verbose && @info "(♻) Predicting time series..."

# take a sub sample of the chain, otherwise the predictions take forever. 
nsamples_chain = min(size(chain, 1), 500)
chain = chain[sample(1:size(chain, 1), nsamples_chain, replace=false), :, :]

# `decondition` model so that `observations` are now treated as random variables.
# This isn't the most elegant way, but it gets the job done

model_decond_ts = DynamicPPL.decondition(model, :obs_ts);
predictions_ts = predict(model_decond_ts, chain)
predictions_ts = [Array(predictions_ts[1:end, :, i]) for i = 1:size(predictions_ts, 3)]

model_decond_ts_prop = DynamicPPL.decondition(model, :obs_ts_prop);
predictions_ts_prop = predict(model_decond_ts_prop, chain)
predictions_ts_prop = [Array(predictions_ts_prop[1:end, :, i]) for i = 1:size(predictions_ts_prop, 3)]

verbose && @info "(✓) Done!"
verbose && @info "(♻) Predicting hospital data..."

model_decond_hosp = DynamicPPL.decondition(model, :obs_hosp);
predictions_hosp = predict(model_decond_hosp, chain)
predictions_hosp = [Array(predictions_hosp[1:end, :, i]) for i = 1:size(predictions_hosp, 3)]

model_decond_hosp_prop = DynamicPPL.decondition(model, :obs_hosp_prop);
predictions_hosp_prop = predict(model_decond_hosp_prop, chain)
predictions_hosp_prop = [Array(predictions_hosp_prop[1:end, :, i]) for i = 1:size(predictions_hosp_prop, 3)]

verbose && @info "(✓) Done!"
verbose && @info "(♻) Predicting seroconversion data..."

model_decond_sero1 = DynamicPPL.decondition(model, :obs_noconv);
predictions_sero1 = predict(model_decond_sero1, chain)
predictions_sero1 = [Array(predictions_sero1[1:end, :, i]) for i = 1:size(predictions_sero1, 3)]

verbose && @info "(✓) Done!"
verbose && @info "(♻) Preparing plots..."

# p2
p2 = plot(; legend=:topleft)
for preds in predictions_ts
    plot_quantiles_ts!(p2, eachcol(preds), 0.95)
end
plot!(p2, data_ts, color=:black, label="observations")
title!(p2, "Weekly incidence")

# p3
p3 = plot(; legend=false)
for preds in predictions_ts_prop
    plot_quantiles!(p3, eachcol(preds), 0.95)
end
bar!(p3, data_ts_prop, color=:red, alpha=0.5)
title!(p3, "Proportional incidence")

# p4
p4 = plot(; legend=false)
for preds in predictions_hosp
    plot_quantiles!(p4, eachcol(preds), 0.95)
end
bar!(p4, data_hosp, color=:red, alpha=0.5)
title!(p4, "Quarterly hospitalizations")

# p5
p5 = plot(; legend=false)
for preds in predictions_hosp_prop
    plot_quantiles!(p5, eachcol(preds), 0.95)
end
bar!(p5, data_hosp_prop, color=:red, alpha=0.5)
title!(p5, "Proportional hospitalizations")

# p6
p6 = plot(; legend=false)
for preds in predictions_sero1
    plot_quantiles!(p6, eachcol(preds), 0.95)
end
bar!(p6, data_noconv, color=:red, alpha=0.5)
title!(p6, "Seroconversion 0-12 months")


# Combine
p_full = plot(
    p_internals, 
    p1, 
    plot(p2, p4, p5, p6, p3, layout= @layout [a ; b; c; d; e]),
    layout=(1,3),
    size=(5000, 3500)
)


@info "Saving plots..."

save(string("output/", filename), p_full)

@info "(✓) Finito!"