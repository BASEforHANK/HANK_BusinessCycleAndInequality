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
    indexes::IndexStruct, compressionIndexes::Array{Array{Int,1},1}, distrSS::Array{Float64,3}; estim=false)
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
    DCD[1]  = mydctmx(n_par.nm_copula)
    DCD[2]  = mydctmx(n_par.nk_copula)
    DCD[3]  = mydctmx(n_par.ny_copula)
    IDCD    = [DCD[1]', DCD[2]', DCD[3]']


    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################
    length_X0 = n_par.ntotal 
    X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,0.0)
    F  = Fsys(X0,X0,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC, IDC,DCD, IDCD)
   
    FR=realpart.(F)
    println("Number of States and Controls")
    println(indexes.profits)
    println("Max error on Fsys:")
    println(maximum(abs.(FR[:])))
    println("Max error of COP in Fsys:")
    println(maximum(abs.(FR[indexes.COP])))
    println("Max error of Vm in Fsys:")
    println(maximum(abs.(FR[indexes.Vm])))
    println("Max error of Vk in Fsys:")
    println(maximum(abs.(FR[indexes.Vk])))
    
    
    ############################################################################
    # Calculate Jacobians of the Difference equation F
    ############################################################################
    BA          = ForwardDiff.jacobian(x-> Fsys(x[1:length_X0],x[length_X0+1:end],XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC, IDC,DCD, IDCD),zeros(2*length_X0))

    B          = BA[:,1:length_X0]
    A          = BA[:,length_X0+1:end]

    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_sgu, nk = SolveDiffEq(A,B, n_par, estim)

    println("State Space Solution Done")

    return gx, hx, alarm_sgu, nk, A, B
end

