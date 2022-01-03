function updateV(EVk::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 c_a_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 c_n_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 m_n_star::Array,#{ForwardDiff.Dual{Nothing,Float64,5},3},
                 r,#::Union{Float64,DualNumbers.Dual{Float64}},
                 q,#::Union{Float64,DualNumbers.Dual{Float64}},
                 m_par::ModelParameters,#::Union{Float64,DualNumbers.Dual{Float64}},
                 n_par::NumericalParameters#::Union{Float64,DualNumbers.Dual{Float64}},
                )
    β::Float64 = m_par.β
    n = size(c_n_star)
    #----------------------------------------------------------------------------
    ## Update Marginal Value Bonds
    #----------------------------------------------------------------------------
    mutil_c_n = mutil(c_n_star)                          # marginal utility at consumption policy no adjustment
    mutil_c_a = mutil(c_a_star)                          # marginal utility at consumption policy adjustment
    Vm   = m_par.λ .* mutil_c_a .+ (1.0 - m_par.λ ) .* mutil_c_n # Expected marginal utility at consumption policy (w &w/o adjustment)
    
    #----------------------------------------------------------------------------
    ## Update marginal Value of Capital
    ## i.e. linear interpolate expected future marginal value of capital using savings policy
    ## Then form expectations.
    #----------------------------------------------------------------------------
    
    Vk = Array{eltype(EVk),3}(undef,n)                          # Initialize Vk-container

@inbounds @views begin
        for j::Int = 1:n[3]
            for k::Int = 1:n[2]
                Vk[:,k,j] = mylinearinterpolate(n_par.grid_m, EVk[:,k,j], m_n_star[:,k,j]) # evaluate marginal value at policy
            end
        end
    end
    Vk        = r .* Vm .+ m_par.λ .* q .* mutil_c_a .+ (1.0 .- m_par.λ) .* β .* Vk # Expected marginal utility at consumption policy (w &w/o adjustment)
    return Vk, Vm
end
