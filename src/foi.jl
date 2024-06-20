# Cosine forcing
function cosine(t, φ, η, β)
    beta_eff = β * (1.0 + η * cos(2.0 *π*(t - φ * 364.0) / 365.0))
    return beta_eff
end
