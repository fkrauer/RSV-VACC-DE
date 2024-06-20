"""
    replace_args(model::DynamicPPL.Model; kwargs...)

Return new `Model` with arguments specified in `kwargs` replaced in `model`.
"""
replace_args(model::DynamicPPL.Model; kwargs...) = replace_args(model, NamedTuple(kwargs))
function replace_args(model::DynamicPPL.Model, args)
    return DynamicPPL.Model{DynamicPPL.getmissings(model)}(
        #=model.name,=# model.f, merge(model.args, args), model.defaults
    )
end

"""
    setmissings(model::DynamicPPL.Model, syms::Symbol...)

Return new `Model` with `syms` set as `missing`.
"""
function setmissings(model::DynamicPPL.Model, syms::Symbol...)
    return DynamicPPL.Model{syms}(model.name, model.f, model.args, model.defaults)
end


"""
    fixedparameters(problem::ODEProblem, prior::Model)

Return `LArray` of parameters from `problem` considered fixed by `prior`.

"""

function fixedparameters(θ, prior::Turing.Model)
    θ_free = prior()
    θ_fixed_labels = ((k for k in keys(θ) if k ∉ keys(θ_free))...,)
    #return @LVector θ[collect(θ_fixed_labels)] θ_fixed_labels
    #vec = @SLVector θ_fixed_labels
    #return vec(θ[collect(θ_fixed_labels)])
    return @LArray θ[collect(θ_fixed_labels)] θ_fixed_labels

end


#=
function fixedparameters(problem::ODEProblem, prior::Turing.Model)
    θ_free = prior()
    θ_fixed_labels = ((k for k in keys(problem.p) if k ∉ keys(θ_free))..., )
    return @LArray problem.p[collect(θ_fixed_labels)] θ_fixed_labels
end
=#

