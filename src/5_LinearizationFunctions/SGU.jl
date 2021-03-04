@doc raw"""
    SGU(XSS,A,B,m_par,n_par,indexes,Copula,compressionIndexes,distrSS;estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys()`](@ref), using Schmitt-Grohé & Uribe (JEDC 2004) style linearization
(apply the implicit function theorem to obtain linear observation and
state transition equations).

The Jacobian is calculated using the package `ForwardDiff`

# Arguments
- `XSS`: steady state around which the system is linearized
- `A`,`B`: matrices to be filled with first derivatives (see `Returns`)
- `m_par::ModelParameters`, `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm
- `Copula::Function`,`distrSS::Array{Float64,3}`: `Copula` maps marginals to
    linearized approximation of joint distribution around `distrSS`
- `indexes::IndexStruct`,`compressionIndexes`: access states and controls by name
    (DCT coefficients of compressed ``V_m`` and ``V_k`` in case of
    `compressionIndexes`)

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_sgu`,`nk`: `alarm_sgu=true` when solving algorithm fails, `nk` number of
    predetermined variables
- `A`,`B`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
"""
function SGU(XSS::Array,A::Array,B::Array, m_par::ModelParameters, n_par::NumericalParameters,
    indexes::IndexStruct, Copula::Function, compressionIndexes::Array{Array{Int,1},1}, distrSS::Array{Float64,3}; estim=false)
    ############################################################################
    # Prepare elements used for uncompression
    ############################################################################
    # Matrices to take care of reduced degree of freedom in marginal distributions
    Γ  = shuffleMatrix(distrSS, n_par)
    # Matrices for discrete cosine transforms
    DC = Array{Array{Float64,2},1}(undef,3)
    DC[1]  = mydctmx(n_par.nm)
    DC[2]  = mydctmx(n_par.nk)
    DC[3]  = mydctmx(n_par.ny)
    IDC    = [DC[1]', DC[2]', DC[3]']

    DCD = Array{Array{Float64,2},1}(undef,3)
    DCD[1]  = mydctmx(n_par.nm-1)
    DCD[2]  = mydctmx(n_par.nk-1)
    DCD[3]  = mydctmx(n_par.ny-1)
    IDCD    = [DCD[1]', DCD[2]', DCD[3]']

    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################
    length_X0 = indexes.profits # Convention is that profits is the last control
    # X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,tuple(zeros(5)...))
    # F  = Fsys(X0,X0,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC, IDC, DCD, IDCD)
    # if maximum(abs.(F))/10>n_par.ϵ
    #     @warn  "F=0 is not at required precision"
    # end
    
    ############################################################################
    # Calculate Jacobians of the Difference equation F
    ############################################################################

    nxB         = length_X0 -length(indexes.Vm) - length(indexes.Vk) 
    nxA         = length_X0 - length(indexes.distr_y) - length(indexes.distr_m) - length(indexes.distr_k)
    BA          = ForwardDiff.jacobian(x-> Fsys(
                                [x[1:indexes.Vm[1]-1]; zeros(length(indexes.Vm) + length(indexes.Vk)); x[indexes.Vm[1]:nxB]],
                                [zeros(length(indexes.distr_y) + length(indexes.distr_m) + length(indexes.distr_k)); x[nxB+1:end]] ,
                                XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC, IDC, DCD, IDCD),zeros(nxB+nxA))
    
    B[:,1:indexes.Vm[1]-1]          = BA[:,1:indexes.Vm[1]-1]
    B[:,indexes.Vk[end]+1:end]      = BA[:,indexes.Vm[1]:nxB]
    A[:,indexes.distr_y[end]+1:end] = BA[:,nxB+1:end]

    # Make use of the fact that Vk/Vm has no influence on any variable in
    # the system, thus derivative is 1
    for i in indexes.Vm
        B[i, i] = 1.0
    end
    for i in indexes.Vk
        B[i, i] = 1.0
    end

    # Make use of the fact that future distribution has no influence on any variable in
    # the system, thus derivative is Γ
    for count = 1:n_par.nm-1 #in eachcol(A)
        i = indexes.distr_m[count]
        A[indexes.distr_m,i] = -Γ[1][1:end-1,count]
    end
    for count = 1:n_par.nk-1 #in eachcol(A)
        i = indexes.distr_k[count]
        A[indexes.distr_k,i] = -Γ[2][1:end-1,count]
    end
    for count = 1:n_par.ny-1 #in eachcol(A)
        i = indexes.distr_y[count]
        A[indexes.distr_y,i] = -Γ[3][1:end-1,count]
    end
    
    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_sgu, nk = SolveDiffEq(A,B, n_par, estim)
    return gx, hx, alarm_sgu, nk, A, B
end

