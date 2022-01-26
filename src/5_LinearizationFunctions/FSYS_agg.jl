@doc raw"""
    Fsys_agg(X, XPrime, Xss, distrSS, m_par, n_par, indexes)

Return deviations from aggregate equilibrium conditions.

`indexes` can be both `IndexStruct` or `IndexStructAggr`; in the latter case
(which is how function is called by [`SGU_estim()`](@ref)), variable-vectors
`X`,`XPrime`, and `Xss` only contain the aggregate variables of the model.
"""
function Fsys_agg(X::AbstractArray, XPrime::AbstractArray, Xss::Array{Float64,1},distrSS::AbstractArray, m_par::ModelParameters,
              n_par::NumericalParameters, indexes::Union{IndexStructAggr,IndexStruct})
              # The function call with Duals takes
              # Reserve space for error terms
    F = zeros(eltype(X),size(X))
    ############################################################################
    #            I. Read out argument values                                   #
    ############################################################################

    ############################################################################
    # I.1. Generate code that reads aggregate states/controls
    #      from steady state deviations. Equations take the form of:
    # r       = exp.(Xss[indexes.rSS] .+ X[indexes.r])
    # rPrime  = exp.(Xss[indexes.rSS] .+ XPrime[indexes.r])
    ############################################################################

    # @generate_equations(aggr_names)
    @generate_equations()

    # Take aggregate model from model file
    # aggregate model marker
    # @include "../3_Model/input_aggregate_model.jl"

    return F
end


