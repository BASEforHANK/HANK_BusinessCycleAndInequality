@doc raw"""
    Kdiff(K_guess,n_par,m_par)

Calculate the difference between the capital stock that is assumed and the capital
stock that prevails under that guessed capital stock's implied prices when
households face idiosyncratic income risk (Aiyagari model).

Requires global functions `employment(K,A,m_par)`, `interest(K,A,N,m_par)`,
`wage(K,A,N,m_par)`, `output(K,TFP,N,m_par)`, and [`Ksupply()`](@ref).

# Arguments
- `K_guess::Float64`: capital stock guess
- `n_par::NumericalParameters`, `m_par::ModelParameters`
"""
function Kdiff(K_guess::Float64, n_par::NumericalParameters, m_par::ModelParameters,
    initial::Bool = true, Vm_guess::AbstractArray = zeros(1, 1, 1),
    Vk_guess::AbstractArray = zeros(1, 1, 1), distr_guess::AbstractArray = zeros(1, 1, 1))

    N       = employment(K_guess, 1.0 ./ (m_par.μ * m_par.μw), m_par)
    r       = interest(K_guess, 1.0 ./ m_par.μ, N, m_par)
    w       = wage(K_guess, 1.0 ./m_par.μ, N, m_par)
    profits = (1.0 .- 1.0 ./ m_par.μ) .* output(K_guess, 1.0, N, m_par)
    RB      = m_par.RB ./ m_par.π
    # K::Float64   = Ksupply(m_par.RB./m_par.π,1.0+r,w*N/n_par.H,profits, n_par,m_par)[1]

    #----------------------------------------------------------------------------
    # Initialize policy function guess
    #----------------------------------------------------------------------------
    # inc[1] = labor income , inc[2] = rental income,
    # inc[3]= liquid assets income, inc[4] = capital liquidation income

    Paux    = n_par.Π^1000
    distr_y = Paux[1, :]
    inc = Array{Array{Float64, 3}}(undef, 4)
    mcw = 1.0 ./ m_par.μw

    # labor income
    # labor income
    incgross = n_par.grid_y .* mcw .* w  .* N ./ n_par.H
    incgross[end] = n_par.grid_y[end] .* profits
    incnet   = m_par.τ_lev .* (mcw .* w .* N ./ n_par.H .* n_par.grid_y).^(1.0 - m_par.τ_prog)
    incnet[end] = m_par.τ_lev .* (n_par.grid_y[end] .* profits).^(1.0 - m_par.τ_prog)
    av_tax_rate = dot((incgross - incnet), distr_y) ./ dot(incgross, distr_y)

    GHHFA = ((m_par.γ - m_par.τ_prog) / (m_par.γ + 1.0)) # transformation (scaling) for composite good
    inc[1] = GHHFA .* m_par.τ_lev .* (n_par.mesh_y .* mcw .* w .* N ./ n_par.H).^(1.0 - m_par.τ_prog) .+
         (1.0 .- mcw) .* w .* N .* (1.0 .- av_tax_rate)# labor income net of taxes
    inc[1][:, : ,end] = m_par.τ_lev .* (n_par.mesh_y[:, :, end] * profits).^(1.0 - m_par.τ_prog) # profit income net of taxes
    # rental income
    inc[2] = r .* n_par.mesh_k
    # liquid asset Income
    eff_int = (RB .+ m_par.Rbar .* (n_par.mesh_m .<= 0.0)) # effective rate
    inc[3] = eff_int .* n_par.mesh_m
    # capital liquidation Income (q=1 in steady state)
    inc[4] = n_par.mesh_k

    # initial guess consumption and marginal values (if not set)
    if initial
        c_guess = inc[1] .+ inc[2].*(n_par.mesh_k .* r .>0) .+ inc[3].*(n_par.mesh_m.>0)
        if any(any(c_guess .< 0.0))
            @warn "negative consumption guess"
        end
        Vm      = eff_int .* mutil(c_guess, m_par.ξ)
        Vk      = (r + m_par.λ) .* mutil(c_guess, m_par.ξ)
        distr   = n_par.dist_guess
    else
        Vm      = Vm_guess
        Vk      = Vk_guess
        distr   = distr_guess
    end

    KS    = Ksupply(RB, 1.0 + r, n_par, m_par, Vm, Vk, distr, inc, eff_int)
    K     = KS[1]
    Vm    = KS[end-2]
    Vk    = KS[end-1]
    distr = KS[end]
    diff    = K - K_guess
    return diff, Vm, Vk, distr
end
