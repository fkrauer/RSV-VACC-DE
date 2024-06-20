function NegativeBinomial2(ψ, incidence; check_args=true)
    p = 1.0/(1.0 + ψ*incidence)
    #p = 1.0/(1.0 + incidence/ψ)
    r = 1.0/ψ
    #r = ψ
    return NegativeBinomial(r, p; check_args=check_args)
end
