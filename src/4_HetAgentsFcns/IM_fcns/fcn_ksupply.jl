@doc raw"""
    Ksupply(RB_guess,R_guess,w_guess,profit_guess,n_par,m_par)

Calculate the aggregate savings when households face idiosyncratic income risk.

Idiosyncratic state is tuple ``(m,k,y)``, where
``m``: liquid assets, ``k``: illiquid assets, ``y``: labor income

# Arguments
- `R_guess`: real interest rate illiquid assets
- `RB_guess`: nominal rate on liquid assets
- `w_guess`: wages
- `profit_guess`: profits
- `n_par::NumericalParameters`
- `m_par::ModelParameters`

# Returns
- `K`,`B`: aggregate saving in illiquid (`K`) and liquid (`B`) assets
-  `TransitionMat`,`TransitionMat_a`,`TransitionMat_n`: `sparse` transition matrices
    (average, with [`a`] or without [`n`] adjustment of illiquid asset)
- `distr`: ergodic steady state of `TransitionMat`
- `c_a_star`,`m_a_star`,`k_a_star`,`c_n_star`,`m_n_star`: optimal policies for
    consumption [`c`], liquid [`m`] and illiquid [`k`] asset, with [`a`] or
    without [`n`] adjustment of illiquid asset
- `V_m`,`V_k`: marginal value functions
"""
function Ksupply(RB_guess::Float64, R_guess::Float64, n_par::NumericalParameters,
    m_par::ModelParameters, Vm::AbstractArray, Vk::AbstractArray, distr_guess::AbstractArray,
    inc::AbstractArray, eff_int::AbstractArray)

    #   initialize distance variables
    dist                = 9999.0
    dist1               = dist
    dist2               = dist

    q                   = 1.0       # price of Capital
    #----------------------------------------------------------------------------
    # Iterate over consumption policies
    #----------------------------------------------------------------------------
    count               = 0
    n                   = size(Vm)
    # containers for policies
    m_n_star            = zeros(n)
    m_a_star            = zeros(n)
    k_a_star            = zeros(n)
    c_a_star            = zeros(n)
    c_n_star            = zeros(n)

    while dist > n_par.ϵ && count < 1000 # Iterate consumption policies until converegence
        count           = count + 1
        # Take expectations for labor income change
        EVk             = reshape(reshape(Vk, (n[1] .* n[2], n[3])) * n_par.Π', (n[1], n[2], n[3]))
        EVm             = reshape((reshape(eff_int, (n[1] .* n[2], n[3])) .*
                                reshape(Vm, (n[1].*n[2], n[3]))) * n_par.Π', (n[1], n[2], n[3]))

        # Policy update step
        c_a_star, m_a_star, k_a_star, c_n_star, m_n_star =
            EGM_policyupdate(EVm, EVk, q, m_par.π, RB_guess, 1.0, inc, n_par, m_par, false)

        # marginal value update step
        Vk_new, Vm_new  = updateV(EVk, c_a_star, c_n_star, m_n_star, R_guess - 1.0, q, m_par, n_par)

        # Calculate distance in updates
        dist1           = maximum(abs, invmutil(Vk_new) .- invmutil(Vk))
        dist2           = maximum(abs, invmutil(Vm_new) .- invmutil(Vm))
        dist            = max(dist1, dist2) # distance of old and new policy

        # update policy guess/marginal values of liquid/illiquid assets
        Vm              = Vm_new
        Vk              = Vk_new
    end
    println(dist)
    
    #------------------------------------------------------
    # Find stationary distribution (Is direct transition better for large model?)
    #------------------------------------------------------

    # Define transition matrix
    S_a, T_a, W_a, S_n, T_n, W_n    = MakeTransition(m_a_star,  m_n_star, k_a_star, n_par.Π, n_par)
    TransitionMat_a                 = sparse(S_a, T_a, W_a, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    TransitionMat_n                 = sparse(S_n, T_n, W_n, n_par.nm * n_par.nk * n_par.ny, n_par.nm * n_par.nk * n_par.ny)
    TransitionMat                   = m_par.λ .* TransitionMat_a .+ (1.0 .- m_par.λ) .* TransitionMat_n

    # Calculate left-hand unit eigenvector
    aux = real.(eigsolve(TransitionMat', 1)[2][1])
    distr = reshape((aux[:]) ./ sum((aux[:])),  (n_par.nm, n_par.nk, n_par.ny))

    #-----------------------------------------------------------------------------
    # Calculate capital stock
    #-----------------------------------------------------------------------------
    K = sum(distr[:] .* n_par.mesh_k[:])
    B = sum(distr[:] .* n_par.mesh_m[:])
    return K, B, TransitionMat, TransitionMat_a, TransitionMat_n, c_a_star, m_a_star, k_a_star, c_n_star, m_n_star, Vm, Vk, distr
end
