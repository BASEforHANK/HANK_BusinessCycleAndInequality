# Elasticities and steepness from target markups for Phillips Curves
η        = μ / (μ - 1.0) # demand elasticity
κ        = η * (m_par.κ / m_par.μ) * (m_par.μ - 1.0) # implied steepnes of phillips curve
ηw       = μw / (μw - 1.0) # demand elasticity wages
κw       = ηw * (m_par.κw / m_par.μw) * (m_par.μw - 1.0) # implied steepnes of wage phillips curve

# Capital Utilization
MPK_SS   = exp(Xss[indexes.rSS]) - 1.0 + m_par.δ_0
δ_1      = MPK_SS
δ_2      = δ_1 .* m_par.δ_s
# Auxiliary variables
Kserv    = K * u # Effective capital
MPKserv  = mc .* Z .* m_par.α .* (Kserv ./ N) .^(m_par.α - 1.0) # marginal product of Capital
depr     = m_par.δ_0 + δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0 # depreciation

Wagesum         = N*w       # Total wages in economy t
WagesumPrime    = NPrime*wPrime # Total wages in economy t+1

############################################################################
#           Error term calculations (i.e. model starts here)          #
############################################################################

#-------- States -----------#
# Error Term on exogeneous States
# Shock processes
F[indexes.Gshock] = log.(GshockPrime) - m_par.ρ_Gshock * log.(Gshock)
F[indexes.Tprogshock] = log.(TprogshockPrime) - m_par.ρ_Pshock * log.(Tprogshock)
F[indexes.Tlevshock] = log.(TlevshockPrime) - m_par.ρ_τshock * log.(Tlevshock)

F[indexes.Rshock] = log.(RshockPrime) - m_par.ρ_Rshock * log.(Rshock)
F[indexes.Sshock] = log.(SshockPrime) - m_par.ρ_Sshock * log.(Sshock)
# Observables only correlated with
F[indexes.τprog_obs]  = log.(τprog_obs) - m_par.α_τ  * (log.(τprog)-Xss[indexes.τprogSS])

F[indexes.A]      = log.(APrime) - m_par.ρ_A * log.(A)    # (unobserved) Private bond return fed-funds spread (produces goods out of nothing if negative)
F[indexes.Z]      = log.(ZPrime) - m_par.ρ_Z * log.(Z)    # TFP
F[indexes.ZI]     = log.(ZIPrime) - m_par.ρ_ZI * log.(ZI) # Invsestment-good Productivity

F[indexes.μ]      = log.(μPrime) - m_par.ρ_μ * log.(μ)  - (1.0 - m_par.ρ_μ) * log.(m_par.μ) # Process for markup target
F[indexes.μw]     = log.(μwPrime) - m_par.ρ_μw * log.(μw)  - (1.0 - m_par.ρ_μw) * log.(m_par.μw) # process for w-amrkup target(NEW)

F[indexes.σ]      = log.(σPrime) - (m_par.ρ_s * log.(σ) + (1.0 - m_par.ρ_s) *
                    m_par.Σ_n * (log(NPrime) - Xss[indexes.NSS]) + log(Sshock)) # Idiosyncratic income risk (contemporaneous reaction to employment)

# Endogeneous States
F[indexes.Ylag] = log(YlagPrime) - log(Y)
F[indexes.Blag] = log(BlagPrime) - log(B)
F[indexes.Zlag] = log(ZlagPrime) - log(Z)
F[indexes.Glag] = log(GlagPrime) - log(G)
F[indexes.Ilag] = log(IlagPrime) - log(I)
F[indexes.wlag] = log(wlagPrime) - log(w)
F[indexes.Tlag] = log(TlagPrime) - log(T)

F[indexes.qlag] = log(qlagPrime) - log(q)
F[indexes.Nlag] = log(NlagPrime) - log(N)
F[indexes.Clag] = log(ClagPrime) - log(C)
F[indexes.πlag] = log(πlagPrime) - log(π)
F[indexes.σlag] = log(σlagPrime) - log(σ)
F[indexes.rlag] = log(rlagPrime) - log(r)
F[indexes.RBlag] = log(RBlagPrime) - log(RB)
F[indexes.av_tax_ratelag] = log(av_tax_ratelagPrime) - log(av_tax_rate)
F[indexes.τproglag] = log(τproglagPrime) - log(τprog)

# Growth rates
F[indexes.Ygrowth] = log(Ygrowth) - log(Y/Ylag)
F[indexes.Tgrowth] = log(Tgrowth) - log(T/Tlag)
F[indexes.Bgrowth] = log(Bgrowth) - log(B/Blag)
F[indexes.Zgrowth] = log(Zgrowth) - log(Z/Zlag)
F[indexes.Ggrowth] = log(Ggrowth) - log(G/Glag)
F[indexes.Igrowth] = log(Igrowth) - log(I/Ilag)
F[indexes.wgrowth] = log(wgrowth) - log(w/wlag)
F[indexes.mcwwgrowth] = log(mcwwgrowth) - log(mcww/mcwwlag)
F[indexes.mcwwlag] = log(mcwwlagPrime) - log(mcww)

F[indexes.qgrowth] = log(qgrowth) - log(q/qlag)
F[indexes.Ngrowth] = log(Ngrowth) - log(N/Nlag)
F[indexes.Cgrowth] = log(Cgrowth) - log(C/Clag)
F[indexes.πgrowth] = log(πgrowth) - log(π/πlag)
F[indexes.σgrowth] = log(σgrowth) - log(σ/σlag)
F[indexes.τproggrowth] = log(τproggrowth) - log(τprog/τproglag)
F[indexes.rgrowth] = log(rgrowth) - log(r/rlag)
F[indexes.RBgrowth] = log(RBgrowth) - log(RB/RBlag)

N_GAP = employment(K, 1.0 ./ (m_par.μ*m_par.μw), m_par)
Y_GAP = output(K,Z,N_GAP, m_par)

#  Taylor rule and interest rates
F[indexes.RB]   = log(RBPrime) - Xss[indexes.RBSS] -
                  ((1 - m_par.ρ_R) * m_par.θ_π).*log(π) -
                  ((1 - m_par.ρ_R) * m_par.θ_Y) .* log(Y/Y_GAP) -
                  m_par.ρ_R * (log.(RB) - Xss[indexes.RBSS])  - log(Rshock)# Taylor rule

# Tax rule
F[indexes.τprog]   = log(τprog) - m_par.ρ_P * log(τproglag)  - (1.0 - m_par.ρ_P) *
                  (Xss[indexes.τprogSS]) - (1.0 - m_par.ρ_P) * m_par.γ_YP * log(Y/Y_GAP) -
                  (1.0 - m_par.ρ_P) * m_par.γ_BP * (log(B)- Xss[indexes.BSS])  - log(Tprogshock)

tax_prog_scale   = (m_par.γ + m_par.τ_prog)/((m_par.γ + τprog))
inc              = [ τlev.*((n_par.mesh_y/n_par.H).^tax_prog_scale .*mcw.*w.*N./(Ht)).^(1.0-τprog)] # capital liquidation Income (q=1 in steady state)
inc[1][:,:,end] .= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes

incgross =[ ((n_par.mesh_y/n_par.H).^tax_prog_scale .*mcw.*w.*N./(Ht))] # capital liquidation Income (q=1 in steady state)
            incgross[1][:,:,end].= (n_par.mesh_y[:,:,end] .* profits)

taxrev   = incgross[1] .- inc[1]
incgrossaux = incgross[1]

F[indexes.τlev] = av_tax_rate - (distrSS[:]' * taxrev[:])./(distrSS[:]' * incgrossaux[:]) # Union profits are taxed at average tax rate

F[indexes.T]    = log(T) - log(distrSS[:]' * taxrev[:] + av_tax_rate*((1.0 .- mcw).*w.*N))
F[indexes.av_tax_rate]       = log(av_tax_rate) - m_par.ρ_τ * log(av_tax_ratelag)  - (1.0 - m_par.ρ_τ) *(Xss[indexes.av_tax_rateSS]) -
                    (1.0 - m_par.ρ_τ) * m_par.γ_Yτ * log(Y/Y_GAP) -
                    (1.0 - m_par.ρ_τ) * m_par.γ_Bτ * (log(B)- Xss[indexes.BSS])  - log(Tlevshock)


# --------- Controls ------------
# Deficit rule
F[indexes.π]   = log(BgrowthPrime) + m_par.γ_B * (log(B)- Xss[indexes.BSS])  -
                 m_par.γ_Y * (log(Y/Y_GAP))  -
                 m_par.γ_π * log(π) - log(Gshock)

F[indexes.G] = log(G) - log(BPrime + T - RB/π*B)

# Phillips Curve to determine equilibrium markup, output, factor incomes (!!add irrelevant shifters of beta!!)
F[indexes.mc]   = (log.(π)- Xss[indexes.πSS]) -
                 (m_par.β * ((log.(πPrime) - Xss[indexes.πSS]) .* YPrime ./ Y) + κ *(mc - 1 ./ μ ))

# Wage Phillips Curve (NEW)
F[indexes.mcw]  = (log.(πw)- Xss[indexes.πwSS]) - (κw *(mcw - 1 ./ μw) +
                    m_par.β * ((log.(πwPrime) - Xss[indexes.πwSS]) .* WagesumPrime ./ Wagesum))
# worker's wage= mcw * firm's wage
# Wage Dynamics (NEW)
F[indexes.πw]   = log.(w./wlag) - log.(πw./π) # LOM for real wage (what firms pay)

# Capital utilisation
F[indexes.u]    = MPKserv  -  q * (δ_1 + δ_2 * (u - 1.0)) # Optimality condition for utilization

# Prices
F[indexes.r]    = log.(r) - log.(1 + MPKserv * u - q * depr ) # rate on capital

F[indexes.mcww] =  log.(mcww) - log.(mcw * w)   # wages

F[indexes.w]    = log.(w) - log.(wage(Kserv, Z * mc, N, m_par))  # wages (NEW)

F[indexes.profits] = log.(profits) - log.(Y * (1.0 - mc))     # profits: + price setting profits + investment profits missing, but =0 on the margin


F[indexes.q]    = 1.0 - ZI * q * (1.0 - m_par.ϕ / 2.0 * (Igrowth - 1.0)^2.0 - # price of capital investment adjustment costs
                   m_par.ϕ * (Igrowth - 1.0) * Igrowth)  -
                   m_par.β * ZIPrime * qPrime * m_par.ϕ * (IgrowthPrime - 1.0) * (IgrowthPrime)^2.0

# Aggregate Quantities
F[indexes.I]   = KPrime .-  K .* (1.0 .- depr)  .- ZI .* I .* (1.0 .- m_par.ϕ ./ 2.0 .* (Igrowth -1.0).^2.0) # Capital accumulation equation
F[indexes.N]   = log.(N) - log.(((1.0 - τprog) * τlev * (mcw.*w).^(1.0 - τprog)).^(1.0 / (m_par.γ+τprog)).*Ht)              # labor supply (NEW)
F[indexes.Y]   = log.(Y) - log.(Z .* N .^(1.0 .- m_par.α) .* Kserv .^m_par.α)         # production function
F[indexes.C]   = log.(Y .- G .- I .+ (A .- 1.0) .* RB .* B ./ π .- (δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0).*K) .- log(C) # Resource constraint

# Error Term on prices/aggregate summary vars (logarithmic, controls), here difference to SS value
# averages
F[indexes.K]      = log.(K)     - Xss[indexes.KSS]
F[indexes.B]      = log.(B)     - Xss[indexes.BSS]
F[indexes.BY]     = log.(BY)    -  log.(B/Y)
F[indexes.TY]     = log.(TY)    -  log.(T/Y)
