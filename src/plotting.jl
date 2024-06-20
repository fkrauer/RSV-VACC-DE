# Quantiles plotting function
plot_quantiles(args...; kwargs...) = plot_quantiles!(plot(), args...; kwargs...)
function plot_quantiles!(p::Plots.Plot, xs, q=0.95; kwargs...)
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

plot_quantiles_hosp(args...; kwargs...) = plot_quantiles_hosp!(plot(), args...; kwargs...)
function plot_quantiles_hosp!(p::Plots.Plot, xs, q=0.95; kwargs...)
    Δ = (1 - q) / 2
    quantiles = mapreduce(hcat, xs) do x
        quantile(x, [Δ, 0.5, 1 - Δ])
    end
    scatter!(
        p,
        quantiles[2, :];
        yerr=(quantiles[2, :] - quantiles[1, :], quantiles[3, :] - quantiles[2, :]),
        yaxis=(:log),
        label="$(q * 100)% credible intervals",
        kwargs...
    )
    return p
end
