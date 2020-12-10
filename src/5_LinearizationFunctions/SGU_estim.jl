@doc raw"""
    SGU_estim(XSS,A,B,m_par,n_par,indexes_aggr,distrSS;estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys`](@ref), while only differentiating with respect to the
aggregate part of the model, [`Fsys_agg()`](@ref).

The partials of the Jacobian belonging to the heterogeneous agent part of the model
are taken from the full-model derivatives provided as arguments, `A` and `B` (computed
by [`SGU()`](@ref)).

# Arguments
- `XSS`: steady state around which the system is linearized
- `A`,`B`: derivative of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
- `m_par::ModelParameters`, `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm
- `indexes::IndexStruct`,`indexes_aggr::IndexStructAggr`: access aggregate states and controls by name
- `distrSS::Array{Float64,3}`: steady state joint distribution

# Returns
as in [`SGU()`](@ref)
"""
function SGU_estim(XSSaggr::Array, A::Array, B::Array,
    m_par::ModelParameters, n_par::NumericalParameters, indexes::IndexStruct,
    indexes_aggr::IndexStructAggr, distrSS::AbstractArray; estim=true)

    ############################################################################
    # Calculate dericatives of non-lineear difference equation
    ############################################################################

    length_X0   = length(XSSaggr) 
    BA          = ForwardDiff.jacobian(x-> Fsys_agg(x[1:length_X0],x[length_X0+1:end],XSSaggr,distrSS,m_par,n_par,indexes_aggr),zeros(2*length_X0))
    Aa          = BA[:,length_X0+1:end]
    Ba          = BA[:,1:length_X0]

    for k = 1:length(aggr_names)
        if !(any(distr_names.==aggr_names[k]))
            j = getfield(indexes, Symbol(aggr_names[k]))
            for h = 1:length(aggr_names)
                if !(any(distr_names.==aggr_names[h]))
                    i = getfield(indexes, Symbol(aggr_names[h]))
                    A[j,i] = Aa[k,h]
                    B[j,i] = Ba[k,h]
                end
            end
        end
    end


    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_sgu, nk = SolveDiffEq(A,B, n_par, estim)
    return gx, hx, alarm_sgu, nk, A, B
end