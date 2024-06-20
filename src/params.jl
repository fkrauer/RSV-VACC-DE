function load_params(epidata; strategy::Int64=0, verbose=true, sensitivity=false)

    # Vectors in the wrapper
    ϵ = 1.0 ./ (365.0 .* epidata.d_y) # ageing vector
    μ = epidata.deaths_d
    n = epidata.n .* 1.0 # vector of sizes of age groups

    # Matrix of proportions symptomatic
    pₛ = ones(25, 5)
    for i ∈ eachindex(epidata.pa)
        pₛ[i, :] .= 1.0 - epidata.pa[i]
    end
    # assume no asymptomatic infections in high risk infants 0-12 months
    pₛ[1:12, 1] .= 1.0 

    # indices
    iₘ = 20:21 # index of age groups with people in child-bearing age

    # demographics (based on estimated 2019 population numbers)
    b = epidata.b0[1] # average daily births

    # fixed transmission parameters:
    φ = 0.5467033 #phase shift of transmission rate (scaled [0,1]), corresponds to mid-January when temperature is lowest
    σ = 1.0 / 4.0 # duration of incubation period
    α = 0.2 # proportional reduction of infectiousness of asymptomatics relative to symptomatics 
    γ₁ = 9.0 # duration of first infection
    g = 0.74 #0.74 means gamma4 is 3.64 days. 0.563 means gamma4 is 1.6 days.
    δ₂ = 1.0 # proportional reduction of susceptibility at level 2, see Sande 2013
    #pₕ = 0.06295 + 0.0087 # Proportion high risk among all births, term 1 is the average proportion preterms (until including week 35, data from TK for 2018-2021), term 2 is the propo of CHD among term, see https://www.thieme-connect.com/products/ejournals/html/10.1055/s-0030-1254155
    pₕ = 0.026 + 0.0088 # born before SSW 35 + born term with CHD, see https://www.aerzteblatt.de/archiv/131886/Zwei-bis-sechs-Wochen-zu-frueh-geboren-Risiken-fuer-das-weitere-Leben and https://www.thieme-connect.com/products/ejournals/html/10.1055/s-0030-1254155#N10DB7
    # BPD are 0.001% of all live births, but mostly among preterms
    
    # estimated transmission parameters:
    β = 0.21 # baseline (constant) transmission rate
    η = 0.18 # amplitude of transmission rate
    ω = 362.0 # duration of immunity after infection
    ξ = 36.0 # duration of maternal immunity
    z = 0.5 # reduction of duration of maternal immunity in high risk/preterm infants
    δ₃ = 0.99 # proportional reduction of susceptibility at level 3 
    δ₄ = 0.38 # proportional reduction of susceptibility at level 4 
    t₃ = 90.0 # duration of third trimester

    # matrix with susceptiblities
    δ = ones(25, 5)
    δ[:,3] .= δ₂
    δ[:,4] .= δ₃
    δ[:,5] .= δ₄

    # Reporting and hospitalization related parameters

    # fixed
    ψ = 0.02 # dispersion parameter of negative binomial
    qₘₑₐₙ = permutedims(DataFrame(CSV.File(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", "qmean_AGI.csv"))).qmean) # Underreporting of lab data
    agemid = epidata.agemid[1:12]
    hₘₑₐₙ = epidata.hmean[13:end]  # hospitalization probability 1-99 year olds
    κₕ = 3.0 # increased risk of hosp of preterm vs. term

    # estimated: these values are some approximations, we don't know the true value before the fit
    ρ₁ = 0.72 # scaling of qₘₑₐₙ
    ρ₂ = 0.68 # scaling of hhospitalization probability h in 1-99 year olds
    a₁ = 0.1 # lower bound coefficient for exponential decay model of hospitalization probability in 0-12 month olds
    a₂ = 2.0  # upper bound coefficient for exponential decay model of hospitalization probability in 0-12 month olds
    a₃ = 3.9 # rate coefficient for exponential decay model of hospitalization probability in 0-12 month olds

    # matrix of hospitalisation probabilities
    h = zeros(25, 5) 
    hh = exponential.(a₁, a₂, a₃, agemid) # high risk months 0-12
    # upper values cannot exceed 1:
    hh = min.(eltype(hh)(1.0), hh)
    hl_ = hh .* 1.0 / κₕ # low risk months 0-12
    hl = vcat(hl_, hₘₑₐₙ .* ρ₂)
    h[:, :] .= hl
    h[1:12,1] .= hh # increased risk only in the first year of life

    # Load Strategies and Intervention related parameters
    if sensitivity 
        interventions = DataFrame(XLSX.readtable(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", "params_intervention_sensitivity_analysis.xlsx"), "Sheet1"))
    else
        interventions = DataFrame(XLSX.readtable(joinpath(dirname(pathof(RSVVACCDE)), "..", "data", "params_intervention.xlsx"), "Sheet1"))
    end
    @assert strategy in eval(Meta.parse.(names(interventions)[3:end])) "chosen strategy no. does not exist in input file"
    iₚₗ, iₚₐ, pₚ, interval_p, ωₚ, nₚ, PEᵢ, PEₛ, PEₕ, tₚ, dₚ, iᵥₗ, iᵥₐ, tiᵥ, v, ωᵥ, VEᵢ, VEₛ, VEₕ, tᵥ, dᵥ, ωᵥₘ, VEᵢₘ, VEₛₘ, VEₕₘ, vₙ = interventions[2:end, string(strategy)]
    
    # Input plausibility checks
    @assert pₚ <= 1.0 && pₚ >= 0.0 "pₚ (mAB uptake) must be [0,1]"
    @assert interval_p == 30.0 || interval_p == 150.0 "The interval of a single mAB treatment must be either 30 days (Palivzumab) or 150 days (Nirsevimab)"
    @assert ωₚ >= 0.0 "ωₚ (duration of protection) must be >=0"
    @assert tₚ >= 0.0 "tₚ (starting day of the immunisation season) must be >=0"
    @assert nₚ <= 5.0 && nₚ >= 1.0 "nₚ (maximum number of injections per season per infant) must be [1,5]"
    @assert dₚ <= 300.0 && dₚ >= 1.0 "duration of seasonal mAB administration must be between 1.0 and 300.0 days"

    @assert v <= 1.0 && v >= 0.0 "v (vaccination update) must be [0,1]"
    @assert ωᵥ >= 0.0 "ωᵥ (duration of protection by vaccination) must be >=0.0"
    @assert ωᵥₘ >= 0.0 "ωᵥₘ (duration of protection by maternal vaccination) must be >= 0.0"
    @assert tᵥ >= 0.0 "tᵥ (timing of vaccination) must be >= 0.0"
    @assert (dᵥ <= 180.0 && dᵥ >= 0.0 && v > 0.0) || (dᵥ == 0.0 && v == 0.0) "dᵥ (duration of seasonal vaccination window) must be between <= 180.0 for seasonal vaccination **or** 0.0 if no vaccination or year-round vaccination"

    @assert vₙ in ["none", "paediatric", "maternal", "older_adult"] "Vaccine name incorrectly specified, must be one of ['none', 'paediatric', 'maternal', 'older_adult']"

    @assert PEᵢ <= PEₛ && PEₛ <= PEₕ "PEᵢ must be <= PEₛ and PEₛ must be <= PEₕ"
    @assert VEᵢ <= VEₛ && VEₛ <= VEₕ "VEᵢ must be <= VEₛ and VEₛ must be <= VEₕ "
    @assert VEᵢₘ <= VEₛₘ && VEₛₘ <= VEₕₘ "VEᵢₘ must be <= VEₛₘ and VEₛₘ must be <= VEₕₘ"

    (dᵥ == 0.0 && vₙ != "none") ? vₛ = false : vₛ = true

    verbose && @info "The immunisation parameters are currently fixed to the following values:"
    verbose && @info "pₚ = $pₚ, ωₚ = $ωₚ, PEᵢ = $PEᵢ, PEₛ = $PEₛ, PEₕ = $PEₕ, v = $v, ωᵥ = $ωᵥ, VEᵢ = $VEᵢ, VEₛ = $VEₛ, VEₕ = $VEₕ, ωᵥₘ = $ωᵥₘ, VEᵢₘ = $VEᵢₘ, VEₛₘ = $VEₛₘ, VEₕₘ = $VEₕₘ"

    # parameter vector
    θ = SLVector(ψ = ψ, 
                ρ₁ = ρ₁, 
                ρ₂ = ρ₂,
                β = β,
                η = η,
                φ = φ,
                α = α,
                ω = ω,
                ξ = ξ,
                z = z,
                σ = σ,
                b = b,
                γ₁ = γ₁,
                g = g,
                δ₂ = δ₂,
                δ₃ = δ₃,
                δ₄ = δ₄,
                t₃ = t₃,
                pₕ = pₕ,
                a₁ = a₁,
                a₂ = a₂,
                a₃ =a₃,
                κₕ = κₕ, 
                pₚ = pₚ,
                ωₚ = ωₚ,
                interval_p = interval_p, 
                PEᵢ = PEᵢ, 
                PEₛ = PEₛ, 
                PEₕ = PEₕ,
                tₚ = tₚ,
                nₚ = nₚ,
                dₚ = dₚ,
                tiᵥ = tiᵥ,
                v = v,
                ωᵥ = ωᵥ,
                VEᵢ = VEᵢ,
                VEₛ = VEₛ,
                VEₕ = VEₕ,
                VEᵢₘ = VEᵢₘ,
                VEₛₘ = VEₛₘ,
                VEₕₘ = VEₕₘ,
                ωᵥₘ = ωᵥₘ,
                tᵥ = tᵥ,
                dᵥ = dᵥ
                )

    return ϵ, μ, n, pₛ, δ, qₘₑₐₙ, hₘₑₐₙ, h, agemid, iₘ, eval(Meta.parse(iₚₐ)), eval(Meta.parse(iₚₗ)), eval(Meta.parse(iᵥₐ)), eval(Meta.parse(iᵥₗ)), vₙ, θ, vₛ

end