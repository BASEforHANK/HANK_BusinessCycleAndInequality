function util(c::AbstractArray, m_par::ModelParameters)
    if m_par.γ == 1.0
        util = log.(c)
    else
        util = c.^(1.0.-m_par.γ) ./ (1.0.-m_par.γ)
    end
    return util
end
