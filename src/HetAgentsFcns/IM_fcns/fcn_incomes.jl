function incomes(n_par,m_par,KSS,distrSS)

    NSS       = employment(KSS, 1.0 ./ (m_par.μ * m_par.μw), m_par)
    rSS       = interest(KSS, 1.0 / m_par.μ, NSS, m_par)
    wSS       = wage(KSS, 1.0 / m_par.μ, NSS, m_par)
    YSS       = output(KSS, 1.0, NSS, m_par)
    ProfitsSS = (1.0 -1.0 / m_par.μ) .* YSS
    ISS       = m_par.δ_0 * KSS
    RBSS      = m_par.RB
    
    
    eff_int   = (RBSS .+ (m_par.Rbar .* (n_par.mesh_m.<=0.0))) # effective rate

    incgross = Array{Array{Float64, 3}}(undef, 6)
    inc = Array{Array{Float64, 3}}(undef, 6)

    mcw = 1.0 ./ m_par.μw

    incgross =[  (n_par.mesh_y .* (1. /m_par.μw).*wSS.*NSS./n_par.H).+
        (1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS,# labor income (NEW)
        interest(KSS,1.0 / m_par.μ, NSS, m_par).* n_par.mesh_k, # rental income
        eff_int .* n_par.mesh_m, # liquid asset Income
        n_par.mesh_k,
        (1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./n_par.H] # capital liquidation Income (q=1 in steady state)
    incgross[1][:,:,end].= n_par.mesh_y[:,:,end] .* ProfitsSS  # profit income net of taxes
    incgross[5][:,:,end].= n_par.mesh_y[:,:,end] .* ProfitsSS  # profit income net of taxes

    inc =[  (((m_par.γ - m_par.τ_prog)/(m_par.γ+1)).*m_par.τ_lev.*(n_par.mesh_y.*1.0 ./m_par.μw.*wSS.*NSS./n_par.H).^(1.0-m_par.τ_prog)).+
        ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS),# labor income (NEW)
        interest(KSS,1.0 / m_par.μ, NSS, m_par).* n_par.mesh_k, # rental income
        eff_int .* n_par.mesh_m, # liquid asset Income
        n_par.mesh_k,
        m_par.τ_lev.*((1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./n_par.H).^(1.0-m_par.τ_prog).*((1.0 - m_par.τ_prog)/(m_par.γ+1)),
        m_par.τ_lev.*((1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./n_par.H).^(1.0-m_par.τ_prog)] # capital liquidation Income (q=1 in steady state)
    inc[1][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)  # profit income net of taxes
    inc[5][:,:,end].= 0.0
    inc[6][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)  # profit income net of taxes

    taxrev        = incgross[5]-inc[6]
    incgrossaux   = incgross[5]
    av_tax_rateSS = (distrSS[:]' * taxrev[:])./(distrSS[:]' * incgrossaux[:])

# apply taxes to union profits
    inc[1] =(((m_par.γ - m_par.τ_prog)/(m_par.γ+1)).*m_par.τ_lev.*(n_par.mesh_y.*1.0 ./m_par.μw.*wSS.*NSS./n_par.H).^(1.0-m_par.τ_prog)).+ # labor income
        ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS).*(1.0 .- av_tax_rateSS)# labor union income
    inc[1][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)
    inc[6] =(m_par.τ_lev.*(n_par.mesh_y.*1.0 ./m_par.μw.*wSS.*NSS./n_par.H).^(1.0-m_par.τ_prog)).+ # labor income
        ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS).*(1.0 .- av_tax_rateSS)# labor union income
    inc[6][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)


return incgross, inc, NSS, rSS, wSS, YSS, ProfitsSS, ISS, RBSS, taxrev, av_tax_rateSS, eff_int
end