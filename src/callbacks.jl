"""
    ForEach

Defines a callable which calls `foreach` on a collection of callables `fs`.

Example usage when used with AbstractMCMC.jl:
```julia
AbstractMCMC.sample(model, sampler, N; callback=ForEach(callback1, callback2))
```
"""
struct ForEach{Fs}
    fs::Fs
end
ForEach(fs::Tuple) = ForEach{typeof(fs)}(fs)
ForEach(fs::AbstractArray) = ForEach{typeof(fs)}(fs)
ForEach(fs...) = ForEach(fs)

function (fe::ForEach)(args...; kwargs...)
    foreach(fe.fs) do f
        f(args...; kwargs...)
    end
end

"""
    StateHistoryCallback

Defines a callable which simply pushes the `state` onto the `states` container.!

Example usage when used with AbstractMCMC.jl:
```julia
# 1. Create empty container for state-history.
state_history = []
# 2. Sample.
AbstractMCMC.sample(model, sampler, N; callback=StateHistoryCallback(state_history))
# 3. Inspect states.
state_history
```
"""
struct StateHistoryCallback{A}
    states::A
end
StateHistoryCallback() = StateHistoryCallback(Any[])

function (cb::StateHistoryCallback)(rng, model, sampler, sample, state, i; kwargs...)
    push!(cb.states, state)
    return nothing
end
