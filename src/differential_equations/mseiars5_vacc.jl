# ODE model
Base.@kwdef struct MSEIARS5_vacc!{C,E,M,S,H,IM,IPA,IPL,IVA,IVL,B,VN,VS}
    c::C=c # contact matrix
    ϵ::E=ϵ # ageing vector
    μ::M=μ # vector of  deaths in each age group
    pₛ::S=pₛ # matrix of proportion symptomatic by age group and level
    h::H=h # matrix of probability of hospitalization when infected by age group and level
    iₘ::IM=iₘ # index of age groups in R compartment that are in child-bearing age
    iₚₐ::IPA=iₚₐ # index of age groups that receive mAB
    iₚₗ::IPL=iₚₗ # index of levels that receive mAB
    iᵥₐ::IVA=iᵥₐ # index of age groups that receive the vaccine
    iᵥₗ::IVL=iᵥₗ # index of levels that receive the vaccine
    betafunc::B=cosine # function for periodic forcing of FOI
    vₙ::VN=vₙ # which vaccine 
    vₛ::VS=vₛ # seasonal vaccine (yes/no)
end

function (wrapper::MSEIARS5_vacc!)(du, u, p, t)

    # Get parameters in vectors and matrices
    c = wrapper.c
    ϵ = wrapper.ϵ
    μ = wrapper.μ
    pₛ = wrapper.pₛ
    h = wrapper.h
    iₘ = wrapper.iₘ
    iₚₐ = wrapper.iₚₐ
    iₚₗ = wrapper.iₚₗ 
    iᵥₐ = wrapper.iᵥₐ
    iᵥₗ = wrapper.iᵥₗ
    vₙ = wrapper.vₙ 
    vₛ = wrapper.vₛ

    # Create views of current states
    N = @view u[:, :, 1:9, :] # All states with demographic changes
    M₁ = @view u[:, :, 1, :]
    M₂ = @view u[:, :, 2, :]
    S = @view u[:, :, 3, :]
    E = @view u[:, :, 4, :]
    I = @view u[:, :, 5, :]
    A = @view u[:, :, 6, :]
    R₁ = @view u[:, :, 7, :]
    R₂ = @view u[:, :, 8, :]
    W = @view u[:, :, 9, :]

    XP = u[1, 1, 10, 2] # seasonal switch for passive immunisation (1=on, 0=off)
    XV = u[1, 1, 10, 3] # seasonal switch for vaccination (1=on, 0=off)   

    # Get parameters in theta
    @unpack α, β, η, φ, σ, b, γ₁, g, ξ, z, ω, δ₂, δ₃, δ₄, t₃, pₕ, pₚ, PEᵢ, PEₛ, PEₕ, ωₚ, v, ωᵥ, ωᵥₘ, VEᵢ, VEₛ, VEₕ, VEᵢₘ, VEₛₘ, VEₕₘ = p[1]

    # Is passive immunisation switched on? Then the proportion mAB treated is >0.0, otherwise it is 0.0
    XP == 1.0 ? pₚ=pₚ : pₚ=0.0

    # Is vaccination switched on? Then the target coverage is >0.0, otherwise it is 0.0
    XV == 1.0 ? v=v : v=0.0

    # Calculate rates from durations. 
    γ = [1.0 / (γ₁ * g^x) for x in [0,0,1,2,3]]
    #Numerator=2 for Erlang distribution with shape 2:
    ω_ = 2.0 / ω # immunity after infection ^-1
    ωₚ_ = 2.0 / ωₚ # immunity after mAB ^-1
    vₙ == "none" ? ωᵥ_ = 0.0 : ωᵥ_ = 2.0 / ωᵥ # immunity after vaccination ^-1
    vₙ == "maternal" ? ωᵥₘ_ = (2.0 / max((ωᵥₘ - ξ), ξ)) : ωᵥₘ_ = 0.0 # immunity after maternal immunisaton ^-1
    t₃_ = 1.0 / t₃ # third trimester ^-1
    ξₕ_ = 2.0 / (ξ * z) # duration of maternal protection in high risk infants ^-1
    ξ_ = 2.0 / ξ # duration of maternal protection in low risk infants ^-1

    # get caches
    δ = p[2]
    n = p[3]
    n = get_tmp(n, du)
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
    prop_preg = p[22]
    prop_preg = get_tmp(prop_preg, du)
    T3_in = p[23]
    T3_in = get_tmp(T3_in, du)
    T3_out = p[24]
    T3_out = get_tmp(T3_out, du)
    hospitalization = p[25]
    hospitalization = get_tmp(hospitalization, du)

    mAB = p[26]
    mAB = get_tmp(mAB, du)
    W_to_E = p[27]
    W_to_E = get_tmp(W_to_E, du)

    TX_waning1 = p[28]
    TX_waning1 = get_tmp(TX_waning1, du)
    TX_waning2 = p[29]
    TX_waning2 = get_tmp(TX_waning2, du)
    TX_waning2_ = p[30]
    TX_waning2_ = get_tmp(TX_waning2_, du)

    vacc_S_OA = p[31]
    vacc_S_OA = get_tmp(vacc_S_OA, du)
    vacc_S_PED = p[32]
    vacc_S_PED = get_tmp(vacc_S_PED, du)
    vacc_S_MAT = p[33]
    vacc_S_MAT = get_tmp(vacc_S_MAT, du)
    cumul_P = p[34]
    cumul_P = get_tmp(cumul_P, du)
    cumul_V = p[35]
    cumul_V = get_tmp(cumul_V, du)

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

    # proportion of women of child-bearing age who are in the untreated R compartments in the third trimester --> maternal protection
    num_pm = sum(R₁[iₘ, :, 5])[1] + sum(R₂[iₘ, :, 5])[1]
    denom_pm = sum(N[iₘ, :, :, 5:6])[1]

    if denom_pm <= 0.0 || num_pm < 0.0 # due to numerical instability around 0 when summing
        pₘ = 0.0
    else 
        pₘ = num_pm / denom_pm
        @assert pₘ <= 1.0 && pₘ >= 0.0 "proportion maternally protected babies from infection is outside range [0,1]"
    end

    # proportion of women of child-bearing age who are in any vaccinated compartments in the third trimester --> maternal protection
    if vₙ == "maternal" #2
        num_pmv = sum(N[iₘ, :, :, 6])[1]
        if denom_pm <= 0.0 || num_pmv < 0.0 # due to numerical instability around 0 when summing
            pₘᵥ = 0.0
        else 
            pₘᵥ = num_pmv / denom_pm
            @assert pₘᵥ <= 1.0 && pₘᵥ >= 0.0 "proportion maternally protected babies from vaccination is outside range [0,1]"
            @assert pₘᵥ + pₘ <= 1.0  "proportion maternally protected babies (infection and vaccination) is larger than 1"
        end
    else
        pₘᵥ = 0.0
    end

    ### normal risk births (level 1)
    births[1, 2, 1, 1] += (b *  (1.0 - pₕ) * pₘ) # births into M1, stratum 1 = protected due to maternal infection
    births[1, 2, 3, 1] += (b * (1.0 - pₕ) * (1.0 - pₘ - pₘᵥ)) # births into S, stratum 1 = not protected
    births[1, 2, 1, 4] +=  (b * (1.0 - pₕ) * pₘᵥ) # births into M1, stratum 4 = protected due to maternal vaccination

    ### high risk births (mainly preterms, level 0)
    births[1, 1, 1, 1] += (b * pₕ * (pₘ + pₘᵥ)) # births into M1, stratum 1 = protected (reduced) due to maternal infection or maternal vaccination
    births[1, 1, 3, 1] += (b * pₕ * (1.0 - pₘ - pₘᵥ)) # births into S, stratum 1 = not protected

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
            S_to_E[i, j, 1] += (λ[i] * δ[i,j] * S[i, j, 1]) # infections in untreated
            S_to_E[i, j, 2] += (λ[i] * δ[i,j] * (1.0 - PEᵢ) * S[i, j, 2])  ## mAB breakthrough 
            S_to_E[i, j, 3] += (λ[i] * δ[i,j] * (1.0 - VEᵢ) * S[i, j, 3]) ## breakthrough infections in vaccinated 
            S_to_E[i, j, 4] += (λ[i] * δ[i,j] * (1.0 - VEᵢₘ) * S[i, j, 4]) ## breakthrough infections in maternally protected 
            S_to_E[i, j, 5] += (λ[i] * δ[i,j] * S[i, j, 5]) ## infections in third trimester women
            S_to_E[i, j, 6] += (λ[i] * δ[i,j] * (1.0 - VEᵢ) * S[i, j, 6]) ## breakthrough infections in third trimester vaccinated women
        end
    end
    # infections in the waning state of mAB or vaccinated
    fill!(W_to_E, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            W_to_E[i, j, 2] += (λ[i] * δ[i,j] * (1.0 - PEᵢ) * W[i, j, 2])  ## mAB breakthrough infections
            W_to_E[i, j, 3] += (λ[i] * δ[i,j] * (1.0 - VEᵢ) * W[i, j, 3]) ## breakthrough infections in vaccinated 
            W_to_E[i, j, 4] += (λ[i] * δ[i,j] * (1.0 - VEᵢₘ) * W[i, j, 4]) ## breakthrough infections in maternally protected 
            W_to_E[i, j, 6] += (λ[i] * δ[i,j] * (1.0 - VEᵢ) * W[i, j, 6]) ## breakthrough infections in third trimester vaccinated women
        end
    end
    
    ## symptomatic infectiousness
    fill!(E_to_I, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            E_to_I[i, j, 1] += (σ * pₛ[i,j] * E[i, j, 1]) # E to I
            E_to_I[i, j, 2] += (σ * pₛ[i,j] * ((1.0 - PEₛ) / (1.0 - PEᵢ)) * E[i, j, 2])  # EP to IP (breakthrough mAB)
            E_to_I[i, j, 3] += (σ * pₛ[i,j] * ((1.0 - VEₛ) / (1.0 - VEᵢ)) * E[i, j, 3])  # EV to IV (breakthrough vaccinated)
            E_to_I[i, j, 4] +=  (σ * pₛ[i,j] * ((1.0 - VEₛₘ) / (1.0 - VEᵢₘ)) * E[i, j, 4]) # EMV to IMV (breakthrough maternally protected)
            E_to_I[i, j, 5] += (σ * pₛ[i,j] * E[i, j, 5]) # E to I third trimester women
            E_to_I[i, j, 6] += (σ * pₛ[i,j] * ((1.0 - VEₛ) / (1.0 - VEᵢ)) * E[i, j, 6])  # EV to IV (breakthrough vaccinated third trimester
        end
    end
    
    ## asymptomatic infectiousness 
    fill!(E_to_A, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            E_to_A[i, j, 1] += (σ * (1.0 - pₛ[i,j]) * E[i, j, 1]) # E to A
            E_to_A[i, j, 2] += (σ * (1.0 - pₛ[i,j] * ((1.0 - PEₛ) / (1.0 - PEᵢ))) * E[i, j, 2])  ## breakthrough mAB
            E_to_A[i, j, 3] += (σ * (1.0 - pₛ[i,j] * ((1.0 - VEₛ) / (1.0 - VEᵢ))) * E[i, j, 3]) ## breakthrough vaccinated
            E_to_A[i, j, 4] += (σ * (1.0 - pₛ[i,j] * ((1.0 - VEₛₘ) / (1.0 - VEᵢₘ))) * E[i, j, 4])  ## breakthrough maternally protected (EMV to AMV)
            E_to_A[i, j, 5] += (σ * (1.0 - pₛ[i,j]) * E[i, j, 5]) # E to A third trimester
            E_to_A[i, j, 6] += (σ * (1.0 - pₛ[i,j] * ((1.0 - VEₛ) / (1.0 - VEᵢ))) * E[i, j, 6]) ## breakthrough vaccinated third trimester
        end
    end
   
    ## clearance of infection

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

    # Pregnancies

    # Proportion pregnant women in the individual states
    fill!(prop_preg, false)
    denom_prop_preg = sum(N[iₘ, :, :, 1])[1]
    @inbounds for k ∈ states
        for j ∈ levels
            for i ∈ agegroups[iₘ]
                prop_preg[i, j, k] += N[i, j, k, 1] / denom_prop_preg
            end
        end
    end

    fill!(T3_in, false)
    @inbounds for k ∈ states[3:end]
        for j ∈ levels
            for i ∈ agegroups[iₘ]
                T3_in[i, j, k] += prop_preg[i, j, k] * b
            end
        end
    end

    fill!(T3_out, false)
    @inbounds for l ∈ strata[5:6]
        for k ∈ states
            for j ∈ levels
                for i ∈ agegroups
                    T3_out[i, j, k, l] += t₃_ * N[i, j, k, l]
                end
            end
        end
    end


    # immunity 

    ## R1 to R2
    fill!(R1_to_R2, false)
    @inbounds for l ∈ strata[[1;5]] #[[1;5;6]]
        for j ∈ levels
            for i ∈ agegroups
                R1_to_R2[i, j, l] += (ω_ * R₁[i, j, l])
            end
        end
    end

    ## waning out of R2
    fill!(waning_out_R2, false)
    @inbounds for l ∈ strata[[1;5]]
        for j ∈ levels
            for i ∈ agegroups
                waning_out_R2[i, j, l] += (ω_ * R₂[i, j, l])
            end
        end
    end

    ## Transitions between levels and strata through waning effects
    fill!(R2_to_S, false)
    @inbounds for l ∈ strata[[1;5]]
        for i ∈ agegroups
            # from level 1,2,3 to levels 2,3,4
            for j ∈ levels[3:end]
                R2_to_S[i, j, l] += waning_out_R2[i, (j-1), l]
            end
            # from level 0 to level 2
            R2_to_S[i, 3, l] += waning_out_R2[i, 1, l]
            # from level 4 to level 4
            R2_to_S[i, end, l] += waning_out_R2[i, end, l] 
        end
    end

    # maternal immunity
    
    ## transition from M1 to M2
    fill!(M1_to_M2, false)
    # low risk
    @inbounds for l ∈ strata[1:4]
        for i ∈ agegroups
            M1_to_M2[i, 2, l] += (ξ_ * M₁[i, 2, l])
        end
    end
    # high risk
    @inbounds for l ∈ strata[1:4]
        for i ∈ agegroups
            M1_to_M2[i, 1, l] += (ξₕ_ * M₁[i, 1, l])
        end
    end

    ## waning out of M2 into S
    fill!(M2_to_S, false)
    # low risk
    @inbounds for l ∈ strata[1:4]
        for i ∈ agegroups
            M2_to_S[i, 2, l] += (ξ_ * M₂[i, 2, l])
        end
    end
    # high risk
    @inbounds for l ∈ strata[1:4]
        for i ∈ agegroups
            M2_to_S[i, 1, l] += (ξₕ_ * M₂[i, 1, l])
        end
    end


    # Prophylactic treatment: mAB
        
    # mAB treatment of those who newly age into the treated age interval
    fill!(mAB, false)
    @inbounds for k ∈ states[1:3]
        for j ∈ levels[iₚₗ]
            mAB[iₚₐ[1], j, k] += (pₚ * N[iₚₐ[1]-1, j, k, 1] * ϵ[iₚₐ[1]-1])
        end
    end

    # Vaccination

    ## Older adults

    ### vaccination of S
    fill!(vacc_S_OA, false)
    if vₙ == "older_adult" 
        if vₛ  # seasonal 
            for j ∈ levels[iᵥₗ]
                # vaccinate those who newly enter the eligible age group through ageing
                vacc_S_OA[iᵥₐ[1], j] += ϵ[iᵥₐ[1]-1] * v * S[iᵥₐ[1]-1, j, 1]
            end
        end
    end

    
    ## paediatric

    ### vaccination of S
    fill!(vacc_S_PED, false)
    if vₙ == "paediatric" 
        for j ∈ levels[iᵥₗ]
            # vaccinate those who newly enter the eligible age group through ageing
            vacc_S_PED[iᵥₐ[1], j] += ϵ[iᵥₐ[1]-1] * v * S[iᵥₐ[1]-1, j, 1]
        end
    end

    ## maternal

    ### vaccination of S
    fill!(vacc_S_MAT, false)
    if vₙ == "maternal" 
        for j ∈ levels[iᵥₗ]
            for i ∈ agegroups[iᵥₐ]
                # vaccinate those who newly enter the third trimester stratum
                vacc_S_MAT[i, j] += T3_in[i, j, 3] * v
            end
        end
    end

    # Waning of immunisations according to Erlang distributions with shape 2

    # from treated S to W
    fill!(TX_waning1, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            TX_waning1[i, j, 2] += (ωₚ_ * S[i, j, 2]) #mABs
            TX_waning1[i, j, 3] += (ωᵥ_ * S[i, j, 3]) #vaccinated
            TX_waning1[i, j, 4] += (ωᵥₘ_ * S[i, j, 4]) #maternally vaccinated
            TX_waning1[i, j, 6] += (ωᵥ_ * S[i, j, 6]) #3rd trimester vaccinated women
        end
    end

    # out of W
    fill!(TX_waning2, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            TX_waning2[i, j, 2] += (ωₚ_ * W[i, j, 2]) #mABs
            TX_waning2[i, j, 3] += (ωᵥ_ * W[i, j, 3]) #vaccinated
            TX_waning2[i, j, 4] += (ωᵥₘ_ * W[i, j, 4]) #maternally vaccinated
            TX_waning2[i, j, 6] += (ωᵥ_ * W[i, j, 6]) #3rd trimester vaccinated women
        end
    end

    # from W back into untreated S
    fill!(TX_waning2_, false)
    # passive immunisation through mAB or maternal vaccination: no gain of immunity (back to same level) and return to statum 1
    @inbounds for l ∈ strata[[2;4]]
        for j ∈ levels
            for i ∈ agegroups
                TX_waning2_[i, j, 1] += TX_waning2[i, j, l] 
            end
        end
    end
    # active immunisation: gain of immunity (to next level) and return to stratum 1
    @inbounds for i ∈ agegroups
        # from level 1,2,3 to levels 2,3,4. Level 0 are not actively immunised
        for j ∈ levels[3:end]
            TX_waning2_[i, j, 1] += TX_waning2[i, (j-1), 3]
        end
        # from level 4 to level 4
        TX_waning2_[i, end, 1] += TX_waning2[i, end, 3] 
    end
    # active immunisation in pregnant women: gain of immunity (to next level) and return to stratum 5
    @inbounds for i ∈ agegroups
        # from level 1,2,3 to levels 2,3,4. Level 0 are not actively immunised
        for j ∈ levels[3:end]
            TX_waning2_[i, j, 5] += TX_waning2[i, (j-1), 6]
        end
        # from level 4 to level 4
        TX_waning2_[i, end, 5] += TX_waning2[i, end, 6] 
    end
   

    # Doses 

    ## cumulative mAB doses given to those who newly enter the lowest eligible age groups (book-keeping)
    fill!(cumul_P, false)
    # Beneficial doses
    @inbounds for k ∈ states[1:3] # M1, M2, S
        for j ∈ levels[iₚₗ]
            for i ∈ agegroups[iₚₐ[1]]
                cumul_P[i, j] += mAB[i, j, k]
            end
        end
    end
    # Wasted doses
    # E, A, R1, R2  (I are not immunised because they are sick)
    @inbounds for k ∈ states[[4;6:end]] 
        for j ∈ levels[iₚₗ]
            for i ∈ agegroups[iₚₐ[1]]
                cumul_P[i, j] += (pₚ * N[iₚₐ[1]-1, j, k, 1] * ϵ[iₚₐ[1]-1])
            end
        end
    end

    ## cumulative vaccine doses given to those who newly enter the vaccinated age groups (book-keeping)
    fill!(cumul_V, false)
    # Older adults
    if vₙ == "older_adult" 
        if vₛ # seasonal 
            @inbounds for j ∈ levels[iᵥₗ]
                for i ∈ agegroups[iᵥₐ[1]]
                    # beneficial doses
                    cumul_V[i, j, 3] += vacc_S_OA[i, j]
                    # wasted doses
                    cumul_V[i, j, 3] += E[i-1, j, 1] * v * ϵ[i-1] 
                    #cumul_V[i, j, 3] += I[i-1, j, 1] * v * ϵ[i-1] # visibly sick are not treated
                    cumul_V[i, j, 3] += A[i-1, j, 1] * v * ϵ[i-1] 
                    cumul_V[i, j, 3] += R₁[i-1, j, 1] * v * ϵ[i-1] 
                    cumul_V[i, j, 3] += R₂[i-1, j, 1] * v * ϵ[i-1] 
                    cumul_V[i, j, 3] += W[i-1, j, 1] * v * ϵ[i-1] 
                end
            end
        end
    end
    
    # Pediatric
    if vₙ == "paediatric" 
        @inbounds for j ∈ levels[iᵥₗ]
            for i ∈ agegroups[iᵥₐ[1]]
                # beneficial doses
                cumul_V[i, j, 3] += vacc_S_PED[i, j]
                # wasted doses
                cumul_V[i, j, 3] += M₁[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += M₂[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += E[i-1, j, 1] * v * ϵ[i-1] 
                #cumul_V[i, j, 3] += I[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += A[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += R₁[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += R₂[i-1, j, 1] * v * ϵ[i-1] 
                cumul_V[i, j, 3] += W[i-1, j, 1] * v * ϵ[i-1] 
            end
        end
    end
    # Maternal
    if vₙ == "maternal" 
        @inbounds for j ∈ levels[iᵥₗ]
            for i ∈ agegroups[iᵥₐ]
                # beneficial doses
                cumul_V[i, j, 6] += vacc_S_MAT[i, j]
                # wasted doses
                cumul_V[i, j, 6] += T3_in[i, j, 4] * v # pregnant women in the E state
                cumul_V[i, j, 6] += T3_in[i, j, 6] * v # pregnant women in the A state
                cumul_V[i, j, 6] += T3_in[i, j, 7] * v # pregnant women in the R1 state
                cumul_V[i, j, 6] += T3_in[i, j, 8] * v # pregnant women in the R2 state
                cumul_V[i, j, 6] += T3_in[i, j, 9] * v # pregnant women in the W state

            end
        end
    end

    # Hospitalizations (modelled as a proportion h of the incident cases)
    fill!(hospitalization, false)
    @inbounds for j ∈ levels
        for i ∈ agegroups
            hospitalization[i, j, 1] += E_to_I[i, j, 1] * h[i, j] #  untreated arm
            hospitalization[i, j, 2] += E_to_I[i, j, 2] * h[i, j] * ((1.0 - PEₕ)/(1.0 - PEₛ))  # mAB arm
            hospitalization[i, j, 3] += E_to_I[i, j, 3] * h[i, j] * ((1.0 - VEₕ)/(1.0 - VEₛ)) ## vaccinated arm
            hospitalization[i, j, 4] += E_to_I[i, j, 4] * h[i, j] *  ((1.0 - VEₕₘ)/(1.0 - VEₛₘ)) ## maternally vaccinated
            hospitalization[i, j, 5] += E_to_I[i, j, 5] * h[i, j] #  third trimester unvaccinated
            hospitalization[i, j, 6] += E_to_I[i, j, 6] * h[i, j] * ((1.0 - VEₕ)/(1.0 - VEₛ)) #third trimester vaccinated
        end
    end

    #----------------------------------------------------------------------------------
    
    # Changes in states

    # N
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
    # into third trimester
    @inbounds for k ∈ states
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, k, 1] -= T3_in[i, j, k]
                du[i, j, k, 5] += T3_in[i, j, k]
            end
        end
    end
    # out of third trimester
    @inbounds for k ∈ states
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, k, 5] -= T3_out[i, j, k, 5]
                du[i, j, k, 1] += T3_out[i, j, k, 5]
                du[i, j, k, 6] -= T3_out[i, j, k, 6]
                du[i, j, k, 3] += T3_out[i, j, k, 6]
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
    @inbounds for j ∈ levels[iₚₗ]
        for i ∈ agegroups[iₚₐ]
            du[i, j, 1, 1] -= mAB[i, j, 1]
            du[i, j, 1, 2] += mAB[i, j, 1]
        end
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
    @inbounds for j ∈ levels[iₚₗ]
        for i ∈ agegroups[iₚₐ]
            du[i, j, 2, 1] -= mAB[i, j, 2]
            du[i, j, 2, 2] += mAB[i, j, 2]
        end
    end

    # S
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 3, l] += R2_to_S[i, j, l]
                du[i, j, 3, l] += M2_to_S[i, j, l]
                du[i, j, 3, l] -= S_to_E[i, j, l]
                du[i, j, 3, l] -= TX_waning1[i, j, l]
                du[i, j, 3, l] += TX_waning2_[i, j, l]
            end
        end
    end
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 3, 1] -= mAB[i, j, 3]
            du[i, j, 3, 2] += mAB[i, j, 3]
            du[i, j, 3, 1] -= vacc_S_OA[i, j]
            du[i, j, 3, 1] -= vacc_S_PED[i, j]
            du[i, j, 3, 3] += vacc_S_OA[i, j]
            du[i, j, 3, 3] += vacc_S_PED[i, j]
            du[i, j, 3, 5] -= vacc_S_MAT[i, j] 
            du[i, j, 3, 6] += vacc_S_MAT[i, j] 
        end
    end

    # E
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 4, l] += S_to_E[i, j, l]
                du[i, j, 4, l] += W_to_E[i, j, l]
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
    @inbounds for l ∈ strata[1:4]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 7, 1] += I_to_R1[i, j, l]
                du[i, j, 7, 1] += A_to_R1[i, j, l]
            end
        end
    end
    @inbounds for l ∈ strata[5:6]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 7, 5] += I_to_R1[i, j, l]
                du[i, j, 7, 5] += A_to_R1[i, j, l]
            end
        end
    end
    @inbounds for l ∈ strata[[1;5]]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 7, l] -= R1_to_R2[i, j, l]
            end
        end
    end

    # R2
    @inbounds for l ∈ strata[[1;5]]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 8, l] += R1_to_R2[i, j, l]
                du[i, j, 8, l] -= waning_out_R2[i, j, l]
            end
        end
    end

    # W
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 9, l] -= W_to_E[i, j, l]
                du[i, j, 9, l] += TX_waning1[i, j, l]
                du[i, j, 9, l] -= TX_waning2[i, j, l]
            end
        end
    end

    # XP (switch)
    du[1, 1, 10, 2] = 0.0

    # XV (switch)
    du[1, 1, 10, 3] = 0.0

    
    # CS: Cumulative symptomatic incidence
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 11, l] = E_to_I[i, j, l]
            end
        end
    end
    

    # CT: Cumulative total incidence
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 12, l] = E_to_A[i, j, l] + E_to_I[i, j, l]
            end
        end
    end
    
    # CH: Cumulative hospitalizations
    @inbounds for l ∈ strata
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 13, l] = hospitalization[i, j, l]
            end
        end
    end

    # D: Cumulative doses
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 14, 2] = cumul_P[i, j]
        end
    end
    @inbounds for l ∈ strata[[3;6]]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 14, l] = cumul_V[i, j, l]
            end
        end
    end

    # TX: treated individuals
    @inbounds for j ∈ levels
        for i ∈ agegroups
            du[i, j, 15, 2] = cumul_P[i, j]
        end
    end
    @inbounds for l ∈ strata[[3;6]]
        for j ∈ levels
            for i ∈ agegroups
                du[i, j, 15, l] = cumul_V[i, j, l]
            end
        end
    end

    return nothing
end
