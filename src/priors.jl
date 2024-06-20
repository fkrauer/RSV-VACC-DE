@model function prior()

    ρ₁ ~ truncated(Normal(1.0, 0.2), 0.5, 1.5)
    ρ₂ ~ truncated(Normal(1.0, 0.2), 0.5, 1.5)
    a₁ ~ Uniform(0.01, 0.2) # Lower bound
    a₂ ~ Uniform(0.2, 3.0) # upper bound
    a₃ ~ Uniform(0.0, 5.0) # rate of decay
    β ~ Beta(2.5, 5.0) 
    η ~ Beta(1.0, 1.0) 
    δ₃ ~ Uniform(0.0, 1.0) #Beta(22.829, 3.171) 
    δ₄ ~ Uniform(0.0, 1.0) #Beta(6.117, 12.882)  
    ω ~ truncated(Gamma(14.0, 16.0), 100.0, 365.0) #truncated(InverseGamma(14.0, 1.0/15.0), 1.0/400.0, 1.0/100.0)
    ξ ~ truncated(Gamma(2.5, 30.0), 10.0, 180.0)

    return @LArray [ρ₁, ρ₂, a₁, a₂, a₃, β, η, δ₃, δ₄, ω, ξ] (:ρ₁, :ρ₂, :a₁, :a₂, :a₃, :β, :η, :δ₃, :δ₄, :ω, :ξ)

end