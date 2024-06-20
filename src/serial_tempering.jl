using AdvancedHMC: Random
using AdvancedHMC: AdvancedHMC
using AdvancedHMC: AbstractMCMC

struct SerialTempering{S,F} <: AbstractMCMC.AbstractSampler
    sampler::S
    schedule::F
end

struct SerialTemperingState{S}
    state::S
    iteration::Int
end

function make_tempered_model(model::AdvancedHMC.DifferentiableDensityModel, β)
    return AdvancedHMC.DifferentiableDensityModel(
        Base.Fix1(*, β) ∘ model.ℓπ, Base.Fix1(.*, β) ∘ model.∂ℓπ∂θ
    )
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG, model, sampler::SerialTempering; kwargs...
)
    # Get the initial temperature.
    β = sampler.schedule(1)
    # Step.
    transition, state = AbstractMCMC.step(
        rng, make_tempered_model(model, β), sampler.sampler; kwargs...
    )

    return transition, SerialTemperingState(state, 1)
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG, model, sampler::SerialTempering, state; kwargs...
)
    # Get the temperature.
    β = sampler.schedule(state.iteration)

    # Step.
    transition, inner_state = AbstractMCMC.step(
        rng, make_tempered_model(model, β), sampler.sampler, state.state; kwargs...
    )

    return transition, SerialTemperingState(inner_state, state.iteration + 1)
end
