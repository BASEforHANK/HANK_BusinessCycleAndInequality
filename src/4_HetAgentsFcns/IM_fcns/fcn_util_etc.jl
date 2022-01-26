#----------------------------------------------------------------------------
# Basic Functions: Utility, marginal utility and its inverse, 
#                  Return on capital, Wages, Employment, Output
#---------------------------------------------------------------------------

function util(c::AbstractArray, m_par::ModelParameters)
    if m_par.γ == 1.0
        util = log.(c)
    else
        util = c.^(1.0.-m_par.γ) ./ (1.0.-m_par.γ)
    end
    return util
end


mutil(c)      = 1.0 ./ ((c.*c).*(c.*c)) # c.^ξ#
invmutil(mu)  = 1.0 ./ (sqrt.(sqrt.(mu))) # mu.^(1.0./ξ)#

# Incomes (K:capital, Z: TFP): Interest rate = MPK.-δ, Wage = MPL, profits = Y-wL-(r+\delta)*K
interest(K::Number, Z::Number,N::Number, m_par::ModelParameters) = Z .* m_par.α .* (K ./ N) .^(m_par.α - 1.0) .- m_par.δ_0
wage(K::Number, Z::Number,N::Number, m_par::ModelParameters)     = Z .* (1 - m_par.α) .* (K./N) .^m_par.α

employment(K::Number, Z::Number, m_par::ModelParameters)         = (Z .* (1.0 - m_par.α) .* (m_par.τ_lev .* (1.0 - m_par.τ_prog)).^(1.0 / (1.0 - m_par.τ_prog)) 
                                                                    .* K .^(m_par.α )).^((1.0 - m_par.τ_prog)./(m_par.γ + m_par.τ_prog + (m_par.α) .* (1 - m_par.τ_prog)))
output(K::Number, Z::Number,N::Number, m_par::ModelParameters)   = Z .* K .^(m_par.α) .* N .^(1 - m_par.α)                                                                    
