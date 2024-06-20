function exponential(L, U, k, age)
    return L + U * exp(-k * age)
end


function gompertz(L, U, k, age_i, age)
    return L + (U-L)*exp(-exp(k*(age - age_i)))
end

