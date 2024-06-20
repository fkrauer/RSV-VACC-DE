# Returns a Dataframe of all incident cases (total and each level) for vaccsim
function ODEwrap_vaccsim(modelfunc,
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
                        strategy, 
                        simyears,
                        n_samples,
                        sensitivity)

    # load params under chosen strategy
    iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, vₙ, θ, vₛ = RSVVACCDE.load_params(epidata; strategy=strategy, sensitivity=sensitivity)[11:end]
    @info "seasonal =  $(vₛ)"
    θ_fixed_labels = ((k for k in keys(θ) if k ∉ estparnames)...,)
    theta_fix = @LArray θ[collect(θ_fixed_labels)] θ_fixed_labels

    # define callback set
    cb_set_new = load_cb_intervention(θ, 
                                    iₚₐ, 
                                    iₚₗ,
                                    iᵥₗ,
                                    iᵥₐ;
                                    vₙ = vₙ,
                                    vₛ = vₛ,
                                    simyears = simyears,
                                    modelfunc = modelfunc,
                                    switch = 10)

    inc_out = DataFrame()
    doses_out = DataFrame()
    immunisations_out = DataFrame()

    for k in ProgressBar(1:n_samples)
        # Get initial states and posterior estimates for given draw
        theta_est = @LArray posterior_sample[k,:] estparnames 
        # remake problem and solve
        θ_new = (@SLVector (keys(theta_est)..., keys(theta_fix)...))(vcat(theta_est, theta_fix))

        # Update matrix for delta (susceptibility)
        δ_new = ones(eltype(θ_new), n_ages, n_levels)
        δ_new[:,3] .= θ_new.δ₂
        δ_new[:,4] .= θ_new.δ₃
        δ_new[:,5] .= θ_new.δ₄

        # Update hospitalisation probabilities (rows: age, cols: levels)
        h = zeros(eltype(θ_new), n_ages, n_levels)
        hh = exponential.(θ_new.a₁, θ_new.a₂, θ_new.a₃, agemid) # high risk months 0-12
        # upper values cannot exceed 1:
        hh = min.(eltype(hh)(1.0), hh)
        hl_ = hh .* 1.0 / θ_new.κₕ # low risk months 0-12
        hvec = vcat(hl_, hₘₑₐₙ .* θ_new.ρ₂)
        h[:,:] .= hvec
        h[1:12,1] .= hh # increased risk only in the first year of life in pre-terms

        # Update problem and solve ODEs
        f_new = MSEIARS5_vacc!(; c, ϵ, μ, pₛ, h, iₘ, iₚₐ, iₚₗ, iᵥₐ, iᵥₗ, betafunc, vₙ, vₛ)
        problem_new = ODEProblem(f_new, u0, tspan, (θ_new, δ_new, cache...))

        sol = solve(problem_new,
                    solver;
                    callback=cb_set_new,
                    solvsettings...)

                    
        # Convert to array
        out = diff(Array(sol), dims=2)
        nt = size(out)[2]
        out = reshape(out, (n_ages, n_levels, 5, n_strata, nt))

        # Count incident cases

        ## Symptomatic incidence
        inc_S = out[:,:,1,:,:]
        # sum over level and stratum
        inc_S = dropdims(sum(inc_S, dims=(2,3)), dims=(2,3))
        inc_S = reshape(inc_S, (n_ages*nt), 1)
        inc_S = DataFrame(inc_S, :auto)
        inc_S[!, :age] = repeat([1:1:n_ages;], outer=nt)
        inc_S[!, :year] = repeat([1:1:nt;], inner=n_ages)
        inc_S[!, :incidence] .= "symptomatic"

        ## Overall incidence
        inc_tot = out[:,:,2,:,:]
        # sum over level and stratum
        inc_tot = dropdims(sum(inc_tot, dims=(2,3)), dims=(2,3))
        inc_tot = reshape(inc_tot, (n_ages*nt), 1)
        inc_tot = DataFrame(inc_tot, :auto)
        inc_tot[!, :age] = repeat([1:1:n_ages;], outer=nt)
        inc_tot[!, :year] = repeat([1:1:nt;], inner=n_ages)
        inc_tot[!, :incidence] .= "total"

        ## Hospital incidence
        inc_H = out[:,:,3,:,:]
        # sum over level and stratum
        inc_H = dropdims(sum(inc_H, dims=(2,3)), dims=(2,3))
        inc_H = reshape(inc_H, (n_ages*nt), 1)
        inc_H = DataFrame(inc_H, :auto)
        inc_H[!, :age] = repeat([1:1:n_ages;], outer=nt)
        inc_H[!, :year] = repeat([1:1:nt;], inner=n_ages)
        inc_H[!, :incidence] .= "hospitalised"

        ## combine to single data frame
        inc = DataFrame()
        append!(inc, inc_tot)
        append!(inc, inc_S)
        append!(inc, inc_H)

        rename!(inc, 1:1 .=> ["n"])
        inc[!,:strategy] .= strategy
        inc[!, :vaccine] .= vₙ
        inc[!,:seasonal] .= vₛ
        inc[!,:replicate] .= k

        append!(inc_out, inc)

        # Count doses

        ## mAB
        doses_mAB = out[:, :, 4, 2, :]
        # sum over level
        doses_mAB = dropdims(sum(doses_mAB, dims=(2)), dims=(2))
        doses_mAB = reshape(doses_mAB, (n_ages*nt), 1)
        doses_mAB = DataFrame(doses_mAB, :auto)
        doses_mAB[!, :age] = repeat([1:1:n_ages;], outer=nt)
        doses_mAB[!, :year] = repeat([1:1:nt;], inner=n_ages)
        doses_mAB[!, :intervention] .= "mAB"

        ## vacc
        doses_vacc = out[:,:,4,[3;6],:]
        # sum over levels and stratum
        doses_vacc = dropdims(sum(doses_vacc, dims=(2,3)), dims=(2,3))
        doses_vacc = reshape(doses_vacc, (n_ages*nt), 1)
        doses_vacc = DataFrame(doses_vacc, :auto)
        doses_vacc[!, :age] = repeat([1:1:n_ages;], outer=nt)
        doses_vacc[!, :year] = repeat([1:1:nt;], inner=n_ages)
        doses_vacc[!, :intervention] .= "vaccine"

        ## combine to single data frame
        doses = DataFrame()
        append!(doses, doses_mAB)
        append!(doses, doses_vacc)

        rename!(doses, 1:1 .=> ["n"])
        doses[!,:strategy] .= strategy
        doses[!, :vaccine] .= vₙ
        doses[!,:seasonal] .= vₛ
        doses[!,:replicate] .= k

        append!(doses_out, doses)

        # count immunisations

        ## mAB
        immunisations_mAB = out[:, :, 5, 2, :]
        # sum over level
        immunisations_mAB = dropdims(sum(immunisations_mAB, dims=(2)), dims=(2))
        immunisations_mAB = reshape(immunisations_mAB, (n_ages * nt), 1)
        immunisations_mAB = DataFrame(immunisations_mAB, :auto)
        immunisations_mAB[!, :age] = repeat([1:1:n_ages;], outer=nt)
        immunisations_mAB[!, :year] = repeat([1:1:nt;], inner=n_ages)
        immunisations_mAB[!, :intervention] .= "mAB"

        ## vacc
        immunisations_vacc = out[:, :, 5, [3; 6], :]
        # sum over levels and stratum
        immunisations_vacc = dropdims(sum(immunisations_vacc, dims=(2, 3)), dims=(2, 3))
        immunisations_vacc = reshape(immunisations_vacc, (n_ages * nt), 1)
        immunisations_vacc = DataFrame(immunisations_vacc, :auto)
        immunisations_vacc[!, :age] = repeat([1:1:n_ages;], outer=nt)
        immunisations_vacc[!, :year] = repeat([1:1:nt;], inner=n_ages)
        immunisations_vacc[!, :intervention] .= "vaccine"

        ## combine to single data frame
        immunisations = DataFrame()
        append!(immunisations, immunisations_mAB)
        append!(immunisations, immunisations_vacc)

        rename!(immunisations, 1:1 .=> ["n"])
        immunisations[!, :strategy] .= strategy
        immunisations[!, :vaccine] .= vₙ
        immunisations[!, :seasonal] .= vₛ
        immunisations[!, :replicate] .= k

        append!(immunisations_out, immunisations)
        
    end

    return inc_out, doses_out, immunisations_out

end

