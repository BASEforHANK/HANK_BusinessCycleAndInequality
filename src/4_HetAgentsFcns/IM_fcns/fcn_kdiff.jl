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
    #----------------------------------------------------------------------------
    # Calculate other prices from capital stock
    #----------------------------------------------------------------------------
    N               = employment(K_guess, 1.0 ./ (m_par.μ * m_par.μw), m_par)   # employment
    w               = wage(K_guess, 1.0 ./m_par.μ, N, m_par)                    # wages
    r               = interest(K_guess, 1.0 ./ m_par.μ, N, m_par)               # Return on illiquid asset
    profits         = (1.0 .- 1.0 ./ m_par.μ) .* output(K_guess, 1.0, N, m_par) # Profit income
    RB              = m_par.RB ./ m_par.π                                       # Real return on liquid assets
    eff_int         = (RB .+ m_par.Rbar .* (n_par.mesh_m .<= 0.0))              # effective rate depending on assets
    GHHFA           = ((m_par.γ + m_par.τ_prog) / (m_par.γ + 1.0))              # transformation (scaling) for composite good
    
    #----------------------------------------------------------------------------
    # Array (inc) to store incomes
    # inc[1] = labor income , inc[2] = rental income,
    # inc[3]= liquid assets income, inc[4] = capital liquidation income
    #----------------------------------------------------------------------------
    Paux            = n_par.Π^1000                                              # Calculate ergodic ince distribution from transitions
    distr_y         = Paux[1, :]                                                # stationary income distribution
    inc             = Array{Array{Float64, 3}}(undef, 4)                        # container for income
    mcw             = 1.0 ./ m_par.μw                                           # wage markup

    # gros (labor) incomes
    incgross        = n_par.grid_y .* mcw .* w  .* N ./ n_par.H                 # gross income workers (wages)
    incgross[end]   = n_par.grid_y[end] .* profits                              # gross income entrepreneurs (profits)
    
    # net (labor) incomes
    incnet          = m_par.τ_lev .* (mcw .* w .* N ./ n_par.H .* n_par.grid_y).^(1.0 - m_par.τ_prog)
    incnet[end]     = m_par.τ_lev .* (n_par.grid_y[end] .* profits).^(1.0 - m_par.τ_prog)
    # average tax rate
    av_tax_rate     = dot((incgross - incnet), distr_y) ./ dot(incgross, distr_y)

    inc[1]          = GHHFA .* m_par.τ_lev .* (n_par.mesh_y .* mcw .* w .* N ./ n_par.H).^(1.0 - m_par.τ_prog) .+
                        (1.0 .- mcw) .* w .* N .* (1.0 .- av_tax_rate) .* n_par.HW         # labor income net of taxes incl. union profits
                        
    inc[1][:,:,end] = m_par.τ_lev .* (n_par.mesh_y[:, :, end] * profits).^(1.0 - m_par.τ_prog) # profit income net of taxes
    
    # incomes out of wealth
    inc[2]          = r .* n_par.mesh_k                                         # rental income
    inc[3]          = eff_int .* n_par.mesh_m                                   # liquid asset income
    inc[4]          = n_par.mesh_k                                              # capital liquidation income (q=1 in steady state)

    #----------------------------------------------------------------------------
    # Initialize policy function (guess/stored values)
    #----------------------------------------------------------------------------

    # initial guess consumption and marginal values (if not set)
    if initial
        c_guess     = inc[1] .+ inc[2].*(n_par.mesh_k .* r .>0) .+ inc[3].*(n_par.mesh_m.>0)
        if any(any(c_guess .< 0.0))
            @warn "negative consumption guess"
        end
        Vm          = eff_int .* mutil(c_guess, m_par.ξ)
        Vk          = (r + m_par.λ) .* mutil(c_guess, m_par.ξ)
        distr       = n_par.dist_guess
    else
        Vm          = Vm_guess
        Vk          = Vk_guess
        distr       = distr_guess
    end
    #----------------------------------------------------------------------------
    # Calculate supply of funds for given prices
    #----------------------------------------------------------------------------
    KS              = Ksupply(RB, 1.0 + r, n_par, m_par, Vm, Vk, distr, inc, eff_int)
    K               = KS[1]                                                     # capital
    Vm              = KS[end-2]                                                 # marginal value of liquid assets
    Vk              = KS[end-1]                                                 # marginal value of illiquid assets
    distr           = KS[end]                                                   # stationary distribution  
    diff            = K - K_guess                                               # excess supply of funds
    return diff, Vm, Vk, distr
end
