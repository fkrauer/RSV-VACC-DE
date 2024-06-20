function load_inits(n, pₕ, n_levels::Int64, n_states::Int64, n_strata::Int64, n_inf::Float64)

    # The u matrix is a 4D element with i=age, j=level, k=state, l=stratum
    n_ages = length(n)
    I0 = vcat(zeros(12), n_inf, zeros(12)) # initially infected 1 year olds
    u0 = zeros(n_ages, n_levels, n_states, n_strata)
    u0[1:12, 1, 3, 1] .= (n[1:12] .* pₕ) # S, <1 year, at level 0
    u0[1:12, 2, 3, 1] .= (n[1:12] .* (1.0 - pₕ)) # S, <1 year, at level 1
    u0[13, 2, 3, 1] = (n[13] - n_inf) # S, 1-year, at level 1
    u0[14:25, 2, 3, 1] .= (n[14:25]) # S, >1 year, at level 1
    u0[:, 2, 4, 1] .= I0 # E, 1 year, at level 1
    
    return u0
end