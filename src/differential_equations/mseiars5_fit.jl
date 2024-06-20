# ODE model
Base.@kwdef struct MSEIARS5_fit!{C,E,M,S,IM,IPA,IPL,B}
    c::C=c # contact matrix
    ϵ::E=ϵ # ageing vector
    μ::M=μ # vector of  deaths in each age group
    pₛ::S=pₛ # matrix of proportion symptomatic by age group and level
    iₘ::IM=iₘ # index of age groups in R compartment that are in child-bearing age
    iₚₐ::IPA=iₚₐ # index of age groups that receive mAB
    iₚₗ::IPL=iₚₗ # index of levels that receive mAB
    betafunc::B=cosine # function for periodic forcing of FOI
end

function (wrapper::MSEIARS5_fit!)(du, u, p, t)

    # Create views of current states
    N = @view u[:, :, 1:8, :] # All states with demographic changes
    M₁ = @view u[:, :, 1, :]
    M₂ = @view u[:, :, 2, :]
    S = @view u[:, :, 3, :]
    E = @view u[:, :, 4, :]
    I = @view u[:, :, 5, :]
    A = @view u[:, :, 6, :]
    R₁ = @view u[:, :, 7, :]
    R₂ = @view u[:, :, 8, :]

    XP = u[1, 1, 9, 2] # seasonal switch for passive immunisation (1=on, 0=off)

    # Get parameters in theta
    @unpack α, β, η, φ, σ, b, γ₁, g, ξ, z, ω, δ₂, δ₃, δ₄, pₕ, pₚ, PEᵢ, PEₛ, ωₚ = p[1]

    # Is passive immunisation switched on? Then the proportion mAB treated is >0.0, otherwise it is 0.0
    XP == 1.0 ? pₚ=pₚ : pₚ=0.0

    # Convert durations to rates
    ξₕ = 2.0 / (ξ * z) # duration of maternal protection in high risk infants ^-1
    ξ = 2.0 / ξ # duration of maternal protection in low risk infants ^-1
    γ = [1.0 / (γ₁ * g^x) for x in [0,0,1,2,3]]
    ω = 2.0 / ω
    ωₚ = 1.0 / ωₚ

    # Get parameters in vectors and matrices
    c = wrapper.c
    ϵ = wrapper.ϵ
    μ = wrapper.μ
    pₛ = wrapper.pₛ
    iₘ = wrapper.iₘ
    iₚₐ = wrapper.iₚₐ
    iₚₗ = wrapper.iₚₗ
    δ = p[2]

    # get caches
    n = p[3]
    n = get_tmp(n, du)
    #c = p[4]
    #c = get_tmp(c, du)
    ageing_out = p[4]
    ageing_out = get_tmp(ageing_out, du)
    ageing_in = p[5]
    ageing_in = get_tmp(ageing_in, du)
    births = p[6]
    births = get_tmp(births, du)
    deaths = p[7]
    deaths = get_tmp(deaths, du)
    I_tot = p[8]
    I_tot = get_tmp(I_tot, du)
    A_tot = p[9]
    A_tot = get_tmp(A_tot, du)
    inf_mat = p[10]
    inf_mat = get_tmp(inf_mat, du)
    λ = p[11]
    λ = get_tmp(λ, du)
    S_to_E = p[12]
    S_to_E = get_tmp(S_to_E, du)
    E_to_I = p[13]
    E_to_I = get_tmp(E_to_I, du)
    E_to_A = p[14]
    E_to_A = get_tmp(E_to_A, du)
    I_to_R1 = p[15]
    I_to_R1 = get_tmp(I_to_R1, du)
    A_to_R1 = p[16]
    A_to_R1 = get_tmp(A_to_R1, du)
    M1_to_M2 = p[17]
    M1_to_M2 = get_tmp(M1_to_M2, du)
    M2_to_S = p[18]
    M2_to_S = get_tmp(M2_to_S, du)
    R1_to_R2 = p[19]
    R1_to_R2 = get_tmp(R1_to_R2, du)
    waning_out_R2 = p[20]
    waning_out_R2 = get_tmp(waning_out_R2, du)
    R2_to_S = p[21]
    R2_to_S = get_tmp(R2_to_S, du)
    mAB = p[22]
    mAB = get_tmp(mAB, du)
    mAB_waning = p[23]
    mAB_waning = get_tmp(mAB_waning, du)

    # define sizes of matrices for rates
    agegroups = axes(N, 1)
    levels = axes(N, 2)
    states = axes(N, 3)
    strata = axes(N, 4)

    # calculate population size by age group
    fill!(n, false)
    @inbounds for l ∈ strata
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups
                    n[i] += N[i, j, k, l]
                end
            end
        end
    end

    #= not necessary, population sizes are stable in this model
    # update contact matrix due to changing pop sizes (make it reciprocal)
    fill!(c, false)
    @inbounds for i ∈ agegroups
        for j ∈ agegroups
            c[i, j] += (c_raw[i, j] * n[i] + c_raw[j, i] * n[j]) / (2.0 * n[i])
        end
    end
    =#

    # Demographic transitions

    ## Ageing out
    fill!(ageing_out, false)
    @inbounds for l ∈ strata
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups
                    ageing_out[i, j, k, l] += (ϵ[i] * N[i, j, k, l])
                end
            end
        end
    end

    ## Ageing in from age groups below
    fill!(ageing_in, false)
    @inbounds for l ∈ strata
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups[2:end]
                    ageing_in[i, j, k, l] += ageing_out[(i-1), j, k, l]
                end
            end
        end
    end

    ## Births
    fill!(births, false)

    # proportion of women of child-bearing age who are in the R compartments --> maternal protection
    num_pm = sum(R₁[iₘ, :, 1])[1] + sum(R₂[iₘ, :, 1])[1]

    if num_pm < 0.0 # due to numerical instability around 0 when summing
        pₘ = 0.0
    else 
        pₘ = num_pm / sum(n[iₘ])[1]
    end

    ### normal risk births
    births[1, 2, 1, 1] += (b * pₘ * (1.0 - pₕ)) # births into M1 compartment at level 1 = protected
    births[1, 2, 3, 1] += (b * (1.0 - pₘ) * (1.0 - pₕ)) # births into S compartment at level 1 = not protected

    ### high risk births (mainly preterms)
    births[1, 1, 1, 1] += (b * pₘ * pₕ) # births into M1 compartment at level 0 = protected
    births[1, 1, 3, 1] += (b * (1.0 - pₘ) * pₕ) # births into S compartment at level 0 = not protected
    
    ## Deaths (calculated from births and transitions in and out, to maintain a stable N per age group)
    fill!(deaths, false)
    @inbounds for l ∈ strata
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups
                    deaths[i, j, k, l] += (μ[i] * (N[i, j, k, l] / n[i]))
                end
            end
        end
    end

    # TRANSMISSION

    ## effective per contact transmission probability
    βeff = wrapper.betafunc(t, φ, η, β)

    ## total number of infectious agents by age group
    fill!(I_tot, false)
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                I_tot[i] += I[i, j, l]
            end
        end
    end

    fill!(A_tot, false)
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                A_tot[i] += A[i, j, l]
            end
        end
    end

    ## total number of infectious contacts by age group
    fill!(inf_mat, false)
    @inbounds for i ∈ agegroups
        for j ∈ agegroups
            inf_mat[i, j] += ((I_tot[j] + A_tot[j] * α) * c[i, j] / n[j])
        end
    end

    ## FOI
    fill!(λ, false)
    @inbounds for j ∈ agegroups
        for i ∈ agegroups
            λ[i] += (βeff * inf_mat[i, j])
        end
    end

    # infections
    fill!(S_to_E, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            S_to_E[i, j, 1] += (λ[i] * δ[i,j] * S[i, j, 1]) # S to E
            S_to_E[i, j, 2] += (λ[i] * δ[i,j] * (1.0 - PEᵢ) * S[i, j, 2])  ## mAB breakthrough infections
        end
    end

    ## symptomatic infectiousness
    fill!(E_to_I, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            E_to_I[i, j, 1] += (σ * pₛ[i,j] * E[i, j, 1]) # E to I
            E_to_I[i, j, 2] += (σ * pₛ[i,j] * ((1.0 - PEₛ)/(1.0 - PEᵢ)) * E[i, j, 2])  # EP to IP (breakthrough mAB)
        end
    end
   
    ## asymptomatic infectiousness 
    fill!(E_to_A, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            E_to_A[i, j, 1] += (σ * (1.0 - pₛ[i,j]) * E[i, j, 1]) # E to A
            E_to_A[i, j, 2] += (σ * (1.0 - pₛ[i,j] * ((1.0 - PEₛ)/(1.0 - PEᵢ))) * E[i, j, 2])  ## breakthrough mAB
        end
    end

    ## clearance symptomatics (I to R1)
    fill!(I_to_R1, false)
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                I_to_R1[i, j, l] += (γ[j] * I[i, j, l]) 
            end
        end
    end

    ## clearance asymptomatics (A to R1)
    fill!(A_to_R1, false)
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                A_to_R1[i, j, l] += (γ[j] * A[i, j, l])
            end
        end
    end

 
    # immunity 

    ## R1 to R2
    fill!(R1_to_R2, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            R1_to_R2[i, j] += (ω * R₁[i, j, 1]) 
        end
    end

    ## waning out of R2
    fill!(waning_out_R2, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            waning_out_R2[i, j] += (ω * R₂[i, j, 1]) 
        end
    end


    ## Transitions between levels: waning from R2 into S
    fill!(R2_to_S, false)
    @inbounds for i ∈ agegroups
        # from level 1,2,3 to levels 2,3,4 
        for j ∈ levels[3:end]
            R2_to_S[i, j] += waning_out_R2[i, (j-1)]
        end
        # from level 0 to level 2
        R2_to_S[i, 3] += waning_out_R2[i, 1]
        # from level 4 to level 4
        R2_to_S[i, end] += waning_out_R2[i, end] 
    end

    
    # Maternal immunity

    ## transition from M1 to M2
    fill!(M1_to_M2, false)
    # low risk
    @inbounds for l ∈ strata
        for i ∈ agegroups
            M1_to_M2[i, 2, l] += (ξ * M₁[i, 2, l])
        end
    end
    # high risk
    @inbounds for l ∈ strata
        for i ∈ agegroups
            M1_to_M2[i, 1, l] += (ξₕ * M₁[i, 1, l])
        end
    end

    ## waning out of M2 into S
    fill!(M2_to_S, false)
    # low risk
    @inbounds for l ∈ strata
        for i ∈ agegroups
            M2_to_S[i, 2, l] += (ξ * M₂[i, 2, l])
        end
    end
    # high risk
    @inbounds for l ∈ strata
        for i ∈ agegroups
            M2_to_S[i, 1, l] += (ξₕ * M₂[i, 1, l])
        end
    end

    
    # Prophylactic treatment: mAB
    
    # mAB treatment of those who newly age into the treated age interval
    # --> only level 0 gets mABs in the fitted model
    fill!(mAB, false)
    @inbounds mAB[iₚₐ[1], 1] += (pₚ * M₁[iₚₐ[1]-1, 1, 1] * ϵ[iₚₐ[1]-1])
    @inbounds mAB[iₚₐ[1], 2] += (pₚ * M₂[iₚₐ[1]-1, 1, 1] * ϵ[iₚₐ[1]-1])
    @inbounds mAB[iₚₐ[1], 3] += (pₚ * S[iₚₐ[1]-1, 1, 1] * ϵ[iₚₐ[1]-1])


    # waning from mAB arm into untreated arm (only level 0 gets mABs in fitted model)
    fill!(mAB_waning, false)
    @inbounds for i ∈ agegroups
        mAB_waning[i] += (ωₚ * S[i, 1, 2])
    end


    #----------------------------------------------------------------------------------
    
    # Changes in states

    # N (all states)
    @inbounds for l ∈ strata
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups
                    du[i, j, k, l] = ageing_in[i, j, k, l]
                    du[i, j, k, l] += births[i, j, k, l]
                    du[i, j, k, l] -= ageing_out[i, j, k, l]
                    du[i, j, k, l] -= deaths[i, j, k, l]
                end
            end
        end
    end

    # M1
    @inbounds for l ∈ strata
        for j ∈ levels[1:2]
            for i ∈ agegroups
                du[i, j, 1, l] -= M1_to_M2[i, j, l]
            end
        end
    end
    # mABs --> only level 0 gets mABs in the fitted model
    @inbounds for i ∈ agegroups
            du[i, 1, 1, 1] -= mAB[i, 1]
            du[i, 1, 1, 2] += mAB[i, 1]
    end

    # M2
    @inbounds for l ∈ strata
        for j ∈ levels[1:2]
            for i ∈ agegroups
                du[i, j, 2, l] += M1_to_M2[i, j, l]
                du[i, j, 2, l] -= M2_to_S[i, j, l]
            end
        end
    end
    # mABs --> only level 0 gets mABs in the fitted model
    @inbounds for i ∈ agegroups
        du[i, 1, 2, 1] -= mAB[i, 2]
        du[i, 1, 2, 2] += mAB[i, 2]
    end


    # S
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 3, l] -= S_to_E[i, j, l]
                du[i, j, 3, l] += M2_to_S[i, j, l]
            end
        end
    end
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 3, 1] += R2_to_S[i, j]
        end
    end
    # mABs --> only level 0 gets mABs in the fitted model
    @inbounds for i ∈ agegroups
        du[i, 1, 3, 1] -= mAB[i, 3]
        du[i, 1, 3, 2] += mAB[i, 3]
        du[i, 1, 3, 2] -= mAB_waning[i]
        du[i, 1, 3, 1] += mAB_waning[i]
    end


    # E
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 4, l] += S_to_E[i, j, l]
                du[i, j, 4, l] -= E_to_I[i, j, l]
                du[i, j, 4, l] -= E_to_A[i, j, l]
            end
        end
    end

    # I
   @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 5, l] += E_to_I[i, j, l]
                du[i, j, 5, l] -= I_to_R1[i, j, l]
            end
        end
    end

    # A
   @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 6, l] += E_to_A[i, j, l]
                du[i, j, 6, l] -= A_to_R1[i, j, l]
            end
        end
    end

    # R1
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 7, 1] += A_to_R1[i, j, l]
                du[i, j, 7, 1] += I_to_R1[i, j, l]
            end
        end
    end
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 7, 1] -= R1_to_R2[i, j]
        end
    end

    # R2
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 8, 1] += R1_to_R2[i, j]
            du[i, j, 8, 1] -= waning_out_R2[i, j]
        end
    end

    # XP (switch)
    du[1, 1, 9, 2] = 0.0

    # CS: Cumulative symptomatic incidence
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 10, l] = E_to_I[i, j, l]
            end
        end
    end

    return nothing
end
