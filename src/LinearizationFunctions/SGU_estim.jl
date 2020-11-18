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
    indexes_aggr::IndexStructAggr, distrSS::AbstractArray; estim=estim)

    ############################################################################
    # Prepare elements used for uncompression
    ############################################################################

    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################

    length_X0 = length(XSSaggr) # Convention is that profits is the last control
    Bd = zeros(length_X0, length_X0)
    Ad = zeros(length_X0, length_X0)

    X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,tuple(zeros(n_FD)...))

    F  = Fsys_agg(X0,X0,XSSaggr,distrSS,m_par,n_par,indexes_aggr)

    @make_deriv_estim n_FD
    BLAS.set_num_threads(1)
    #@timev begin
    prime_loop_estim!(Ad, DerivPrime, length_X0, n_par,n_FD)
    loop_estim!(Bd, Deriv, length_X0, n_par,n_FD)
    #end
    #@timev begin
    for k = 1:length(aggr_names)
        if !(any(distr_names.==aggr_names[k]))
            j = getfield(indexes, Symbol(aggr_names[k]))
            for h = 1:length(aggr_names)
                if !(any(distr_names.==aggr_names[h]))
                    i = getfield(indexes, Symbol(aggr_names[h]))
                    A[j,i] = Ad[k,h]
                    B[j,i] = Bd[k,h]
                end
            end
        end
    end
    #end
    #A[n_par.Asel] = Ad
    #B[n_par.Bsel] = Bd

    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    BLAS.set_num_threads(Threads.nthreads())
    if n_par.sol_algo == :lit # naive linear time iteration
        @views begin
            F1 = A[:, 1:n_par.nstates]
            F2 = A[:, n_par.nstates+1:end]
            F3 = B[:, 1:n_par.nstates]
            F4 = B[:, n_par.nstates+1:end]
            BB = hcat(F1, F4)
            AA = F3 #hcat(F3, zeros(n_par.ntotal,n_par.ncontrols))
            CC = hcat(zeros(n_par.ntotal,n_par.nstates),F2)

            F0 = zeros(n_par.ntotal, n_par.nstates)
            F00 = zeros(n_par.ntotal, n_par.nstates)

            F0[1:n_par.nstates, 1:n_par.nstates]     .= n_par.LOMstate_save
            F0[n_par.nstates+1:end, 1:n_par.nstates] .= n_par.State2Control_save

            diff1 = 1000.0
            i = 0
            Mat = copy(BB)
            #@timev begin
            while abs(diff1)>1e-6 && i<1000
                Mat[:,1:n_par.nstates] = BB[:,1:n_par.nstates] .+ CC * F0
                F00 = Mat \ (-AA)
                #F00 .= (BB .+ CC * F0) \ (-AA)
                diff1 = maximum(abs.(F00[:] .- F0[:]))[1]
                F0 .= F00
                i += 1
            end
            #end
            hx = F0[1:n_par.nstates,1:n_par.nstates]
            gx = F0[n_par.nstates+1:end,1:n_par.nstates]
            #println(i)
            #println("lit")
            nk = n_par.nstates
            alarm_sgu = false
        end
    elseif n_par.sol_algo == :litx # linear time iteration with speed-up
         @views begin
            # Formula X_{t+1} = [X11 X12;X21 X22] = inv(BB + [0 F2]*X_t)*[F3 0]
            # implies X21=0, X22 = 0
            # use of Woodburry-Morrison-Shermann
            # X.1_{t+1} =  (inv(BB)  - inv(BB)*[0 F2]* inv(I+[X.1_t 0]*inv(BB)*[0 F2])[X.1_t 0])*F3inv(BB) ; WMS formula
            # Z = inv(BB)
            # X.1_{t+1} =  Z*F3  - Z * [0 F2] * inv(I+[X.1_t 0]*Z*[0 F2])[X.1_t 0] * Z * F3 ; multiply
            #           =  Z*F3  - Z * [0 F2] * inv(I+[0 X11_t *Z1. * F2; 0 X21_t* Z2. *F2])[X.1_t 0]*Z*F3
            #           =  Z*F3  - Z * [0 F2] * [I  , -(X11_t * Z1. * F2) *inv(I+ X21_t* Z2. *F2); 0 , inv(I+ X21_t* Z2. *F2)])[X.1_t 0]*Z*F3 ; block upper trianglular inverse
            #           =  Z*F3  - [0  , Z * F2 * inv(I+ X21_t* Z2. *F2)] * X.1_t * Z1. *F3 ; multiply
            #           =  Z*F3  - Z*F2* inv(I+ X21_t* Z2. *F2)]* X21_t * Z1. *F3 ; multiply out
            #           =  [(Z1. * F3); (Z2. * F3) ] -  [(Z1. * F2); (Z2. * F2)]* (I - X21_t * inv(I +  (Z2. * F2) * X21_t)]* (Z2. * F2)) * X21_t * (Z1. *F3); WMS formula
            # X11_{t+1} =  Q3 -  Q1 *(I - X21_t*inv(I +  Q1*X21_t) * Q1)*X21_t * Q3; WMS formula
            #           =  Q3 -  Q1 *(X21_t - X21_t*inv(I +  Q1*X21_t) * Q1*X21_t) * Q3;
            # X21_{t+1} =  Q4 -  Q2 *(X21_t - X21_t*inv(I +  Q1*X21_t) * Q1*X21_t) * Q3; WMS formula
            #           =  Q4 -  Q2 *X21_t*inv(I +  Q1*X21_t) * Q3; use inv(I+A)*A = I - inv(I+A)  
            #        Q1 = Z1. * F2
            #        Q2 = Z2. * F2
            #        Q3 = Z1. * F3
            #        Q4 = Z2. * F3
            F1 = A[:, 1:n_par.nstates]
            F2 = A[:, n_par.nstates+1:end]
            F3 = -B[:, 1:n_par.nstates]
            F4 = B[:, n_par.nstates+1:end]
            
            BB = hcat(F1, F4)
            Z  = BB\I # inverse of BB

            X21up = zeros(n_par.ncontrols, n_par.nstates)
            
            X21   = copy(n_par.State2Control_save)

            diff1 = 1000.0
            i = 0
            Q1 = Z[1:n_par.nstates,:]*F2    #[Z11 Z12] * F2
            Q2 = Z[n_par.nstates+1:end,:]*F2#[Z21 Z22] * F2
            Q3 = Z[1:n_par.nstates,:]*F3    #[Z11 Z12]*F3
            Q4 = Z[n_par.nstates+1:end,:]*F3#[Z21 Z22]*F3
            Q1X21 = Q1*X21
            #H  = X21 - X21*((I+Q1X21)\Q1X21)
            H     = X21/(I + Q1X21)
            # within loop only update X21
            #@timev begin
            while diff1>1e-6 && i<1000 
                i += 1
                X21up .= Q4 - Q2*H*Q3#  [Z21 - Q2*H*X21*Z11  Z22 - Q2*H*X21*Z12]*F3
                diff1  = maximum(abs.(X21up .- X21))[1]
                Q1X21 .= Q1*X21up
                X21   .= X21up
                H     .= X21/(I + Q1X21) # use that inv(I+A)*A = I - inv(I+A) 
            end
            #end
            X11 = (I + Q1X21)\Q3 # Q3 - Q1*H*Q3 # [Z11 - Q1 * H*X21*Z11  Z12 - Q1*H*X21*Z12]* F3
            hx = X11#F0[1:n_par.nstates,1:n_par.nstates]
            gx = X21#[n_par.nstates+1:end,1:n_par.nstates]
            #println(i)
            #println("litx")
            nk = n_par.nstates
            alarm_sgu = false
        end
    elseif n_par.sol_algo == :aim # Anderson and Moore method
        F1 = A[:, 1:n_par.nstates]
        F2 = A[:, n_par.nstates+1:end]
        F3 = B[:, 1:n_par.nstates]
        F4 = B[:, n_par.nstates+1:end]
        BB = hcat(F1, F4)
        AA = hcat(F3, zeros(n_par.ntotal,n_par.ncontrols))
        CC = hcat(zeros(n_par.ntotal, n_par.nstates), F2)
        FF = [AA BB CC]

        epsi   = 2.2e-16
        condn  = 1.e-10
        uprbnd = 1.0 + 1.e-6
        F0, rts, ia, nexact, nnumeric, lgroots, aimcode =
            AndersonMooreAlg(copy(FF), n_par.ntotal, 1, 1, condn, uprbnd)
        if aimcode == 1
            hx = F0[1:n_par.nstates, 1:n_par.nstates]
            gx = F0[n_par.nstates+1:end, 1:n_par.nstates]
            alarm_sgu = false
        else
            hx = Array{Float64}(undef, n_par.nstates, n_par.nstates)
            gx = Array{Float64}(undef, n_par.ncontrols, n_par.nstates)
            alarm_sgu = true
        end
        nk = n_par.nstates
    elseif n_par.sol_algo == :schur # (complex) schur decomposition
        alarm_sgu = false
        Schur_decomp, slt, nk, λ = complex_schur(A, -B) # first output is generalized Schur factorization

        # Check for determinacy and existence of solution
        if n_par.nstates != nk
            if estim # return zeros if not unique and determinate
                hx = Array{Float64}(undef, n_par.nstates, n_par.nstates)
                gx = Array{Float64}(undef, n_par.ncontrols, n_par.nstates)
                alarm_sgu = true
                return gx, hx, alarm_sgu, nk, A, B
            else # debug mode/ allow IRFs to be produced for roughly determinate system
                ind    = sortperm(abs.(λ); rev=true)
                slt = zeros(Bool, size(slt))
                slt[ind[1:n_par.nstates]] .= true
                alarm_sgu = true
                @warn "critical eigenvalue moved to:"
                print(λ[ind[n_par.nstates - 5:n_par.nstates+5]])
                print(λ[ind[1]])
                nk = n_par.nstates
            end
        end
        # in-place reordering of eigenvalues for decomposition
        ordschur!(Schur_decomp, slt)

        # view removes allocations
        z21 = view(Schur_decomp.Z, (nk+1):n_par.ntotal, 1:nk)
        z11 = view(Schur_decomp.Z, 1:nk, 1:nk)
        s11 = view(Schur_decomp.S, 1:nk, 1:nk)
        t11 = view(Schur_decomp.T, 1:nk, 1:nk)

        if rank(z11) < nk
            @warn "invertibility condition violated"
            hx = Array{Float64}(undef, n_par.nstates, n_par.nstates)
            gx = Array{Float64}(undef, n_par.ncontrols, n_par.nstates)
            alarm_sgu = true
            return gx, hx, alarm_sgu, nk, A, B
        end
        z11i = z11 \ I # I is the identity matrix -> doesn't allocate an array!
        gx = real(z21 * z11i)
        hx = real(z11 * (s11 \ t11) * z11i)
    else
        error("Solution algorithm not defined!")
    end
    return gx, hx, alarm_sgu, nk, A, B
end

# Calculating Jacobian to XPrime
function prime_loop_estim!(Ad, DerivPrime, length_X0, n_par, n_FD)
    Threads.@threads for i = 1:n_FD:length_X0#in eachcol(A)
        aux = DerivPrime.(i)
        for k = 1:min(n_FD,length_X0-i+1)
            for j = 1:size(Ad,1)
                Ad[j,i+k-1] = aux[j][k]
            end
        end
    end
end

# Calculating Jacobian to X
function loop_estim!(Bd, Deriv, length_X0, n_par, n_FD)
    Threads.@threads for i = 1:n_FD:length_X0#in eachcol(A)
        aux = Deriv.(i)
        for k = 1:min(n_FD,length_X0-i+1)
            for j = 1:size(Bd,1)
                Bd[j,i+k-1] = aux[j][k]
            end
        end
    end
end
