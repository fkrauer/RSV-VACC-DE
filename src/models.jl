@model function turingmodel(prior,
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
    sensealg=ForwardDiffSensitivity())

    issuccess = true

    # Sample prior parameters.
    theta_est = @submodel prior()

    # Update `θ` 
    θ = (@SLVector (keys(theta_est)..., keys(theta_fix)...))(vcat(theta_est, theta_fix))
    
    # Update matrix for delta (susceptibility)
    δ = ones(eltype(θ), 25, 5)
    δ[:,3] .= θ.δ₂
    δ[:,4] .= θ.δ₃
    δ[:,5] .= θ.δ₄
    
    # Update problem and solve ODEs
    problem_new = remake(problem; p=(θ, δ, cache...))

    sol = solve(problem_new,
                solver;
                sensealg=sensealg,
                callback=cb_set,
                solvsettings...)

    # Return early if integration failed
    issuccess &= (sol.retcode === :Success) # SciMLBase.successful_retcode(sol)
    if !issuccess
        Turing.@addlogprob! -Inf
        return nothing
    end

    # Convert solution to array
    sol_array = Array(sol)

    # Compute daily symptomatic incidence by age, summed over all levels
    cumul = sol_array[:,:,10,:,:]
    inc_d = diff(dropdims(sum(cumul, dims=(2,3)), dims=(2,3)), dims=2)
    inc_d = transpose(inc_d)

    # Avoid negative incidence due to numerical instability issue.
    if minimum(inc_d) <= eltype(inc_d)(0.0)
        inc_d = max.(eltype(inc_d)(1e-9), inc_d)
    end

    # 1. calculate weekly, age-specific incidence, for comparison to AGI time series
    inc_w = reshape(inc_d, 7, div(size(inc_d,1),7), 25)
    inc_w = dropdims(sum(inc_w, dims=1), dims=1)

    # bin to age groups in the data
    inc_w_rep = hcat(sum(inc_w[:,1:13], dims=2), # 0-2
                        sum(inc_w[:,14:16], dims=2), # 2-5
                        sum(inc_w[:,17:18], dims=2), # 5-15
                        sum(inc_w[:,19:20], dims=2), # 15-35
                        sum(inc_w[:,21:25], dims=2)) # 35+

    # calculate reported incidence
    inc_w_rep = inc_w_rep .* qₘₑₐₙ .* θ.ρ₁
    
    # Remove cases from time points/age groups not observed
    inc_w_rep[setdiff(eachindex(inc_w_rep), ts_obsindex)] .= 0.0

    # Calculate proportion of age groups in the total incidence
    inc_sum = sum(inc_w_rep, dims=1)
    prop_ts = vec(inc_sum ./ sum(inc_sum, dims=2))
   
    # Sum by time point over all age groups
    inc_w_rep_sum = sum(inc_w_rep, dims=2)
   
    # Subset to observed time points
    inc_w_rep_sum = inc_w_rep_sum[convert.(Int, ((times_obsindex .- 4) ./ 7) .+ 1), :]

    # Reshape time series into 1D matrix
    inc_w_rep = reshape(inc_w_rep, (size(inc_w_rep,2) * size(inc_w_rep,1), 1))
   
    # subset to observed time points
    inc_w_rep = inc_w_rep[ts_obsindex,:]
    
    # 2. Age-stratified hospitalisation incidence by quarter to compare to TK data

    # calculate the daily incidence by age and level
    inc_d_level = diff(cumul[:,:,:,1:1461], dims=4)

    # Avoid negative incidence due to numerical instability issue.
    if minimum(inc_d_level) <= eltype(inc_d_level)(0.0)
        inc_d_level = max.(eltype(inc_d_level)(1e-9), inc_d_level)
    end

    # Construct the matrix h that holds the hospitalisation probabilities (rows: age, cols: levels)
    h = zeros(eltype(θ), 25, 5, 2)
    hh = exponential.(θ.a₁, θ.a₂, θ.a₃, agemid) # high risk months 0-12
    # upper values cannot exceed 1:
    hh = min.(eltype(hh)(1.0), hh)
    hl_ = hh .* 1.0 / θ.κₕ # low risk months 0-12
    hl = vcat(hl_, hₘₑₐₙ .* θ.ρ₂)
    h[:, :, :] .= hl
    h[1:12,1,:] .= hh # increased risk only in the first year of life in pre-terms
    # reduction of h due to Palivizumab
    hp = hh[2:12] .* ((1.0 - θ.PEₕ)/(1.0 - θ.PEₛ))
    h[2:12, 1, 2] .= hp
 
    # Sum over strata (treated and untreated) and calculate the incident hospitalised cases
    hosp_d_level = dropdims(sum(inc_d_level .* h, dims=3), dims=3)
    prop_highrisk = hosp_d_level[:,1:1,1] ./ sum(hosp_d_level[:,:,1], dims=2)

    # Sum over all levels
    hosp_d = dropdims(sum(hosp_d_level, dims=2), dims=2)
    hosp_d = permutedims(hosp_d)

    # Aggregate by quarterly cases
    hosp_q = zeros(eltype(hosp_d), length(dim_q), 25) 
    @inbounds for age in 1:25
        start = 1
        rowcount = 1
        for len in dim_q
            hosp_q[rowcount, age] = sum(hosp_d[start:(start + len - 1), age], dims=1)[1]
            start += len
            rowcount += 1
        end
    end
    
    # Remove cases from time points/age groups not observed
    hosp_q[setdiff(eachindex(hosp_q), hosp_obsindex)] .= 0.0

    # Calculate proportion of age groups in the total hospitalised cases
    hosp_sum = sum(hosp_q, dims=1)
    prop_hosp = vec(hosp_sum ./ sum(hosp_sum, dims=2))

    # keep model predictions for age groups with data
    hosp_q_rep = vec(hosp_q) #reshape(hosp_q, (size(hosp_q,1) * size(hosp_q,2), 1))
    hosp_q_rep = hosp_q_rep[hosp_obsindex, :]
    
    # Sum by time point over all age groups
    hosp_q_rep_sum = sum(hosp_q, dims=2)
    
    # Remove rows with zeros (times not observed)
    hosp_q_rep_sum = hosp_q_rep_sum[findall(hosp_q_rep_sum .!= 0.0),:]

    # 3. calculate proportion never infected to compare with seroprevalence data

    # Calculate the average number of individuals in each age groups 2-12 months at each time point
    N = sum(sol_array[2:12, : ,1:8, :, :], dims=(2,3,4))
    N = dropdims(N, dims=(2,3,4))

    # calculate n in age groups 2:12 in compartments M, S, E, I, A at level 1, at each time point
    n_noconv = sum(sol_array[2:12, 1:2, 1:6, :, :], dims=(2,3,4))
    n_noconv = dropdims(n_noconv, dims=(2,3,4))

    # calculate proportion not yet seroconverted
    prop_noconv = n_noconv ./ N
    # calculate average over all time points
    prop_noconv = mean(prop_noconv, dims=2)
    # correct for small numerical errors if p > 1
    prop_noconv = min.(eltype(prop_noconv).(1.0), prop_noconv)

    # Likelihood function.
    # HACK(t-tfjelde): makes it so that when we compute pointwise log-likelihood,
    # we treat the observations as a collection of independent observations rather
    # than a single multivariate observation.
    
    if __context__ isa DynamicPPL.PointwiseLikelihoodContext
        
        obs_ts = Array{Int}(undef, size(inc_w_rep_sum))
        for i in eachindex(inc_w_rep_sum)
            obs_ts[i] ~ NegativeBinomial2(θ.ψ, inc_w_rep_sum[i])
        end
        obs_ts_prop ~ Multinomial(N_ts, prop_ts)    
        
        obs_hosp = Array{Int}(undef, size(hosp_q_rep_sum))
        for i in eachindex(hosp_q_rep_sum)
            obs_hosp[i] ~ NegativeBinomial2(θ.ψ, hosp_q_rep_sum[i])
        end
        obs_hosp_prop ~ Multinomial(N_hosp, prop_hosp)

        obs_noconv = Array{Int}(undef, size(prop_noconv))
        for i in eachindex(prop_noconv)
            obs_noconv[i] ~ Binomial(N_noconv[i], prop_noconv[i])
        end

    else
        
        obs_ts ~ arraydist(@. NegativeBinomial2(θ.ψ, inc_w_rep_sum))
        obs_ts_prop ~ Multinomial(N_ts, prop_ts)
        obs_hosp ~ arraydist(@. NegativeBinomial2(θ.ψ, hosp_q_rep_sum))
        obs_hosp_prop ~ Multinomial(N_hosp, prop_hosp)
        obs_noconv ~ arraydist(@. Binomial(N_noconv, prop_noconv))

    end

    return (; sol,
            θ=θ,
            hl,
            hh,
            inc_w_rep,
            prop_highrisk,
            hosp_q_rep,
            obs_ts,
            obs_ts_prop,
            obs_hosp,
            obs_hosp_prop,
            obs_noconv,
            ), issuccess

end

