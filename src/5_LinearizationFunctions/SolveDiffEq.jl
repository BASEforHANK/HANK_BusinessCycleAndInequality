@doc raw"""
    SolveDiffEq(A, B, n_par; estim)

Calculate the solution to the linearized difference equations defined
B x_t = A x_{t+1}.

# Arguments
- `A`,`B`: matrices with first derivatives 
- `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm, options are: 
    * `lit`:   Linear time iteration (naive implementation, not recommended)
    * `litx`:  Linear time iteration (improved implementation following Reiter fast if initial guess is good)
    * `schur`: Klein's algorithm (fastest if no good initial guess is available)

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_sgu`,`nk`: `alarm_sgu=true` when solving algorithm fails, `nk` number of
    predetermined variables
"""
function SolveDiffEq(A::Array,B::Array, n_par::NumericalParameters, estim=false)
    if n_par.sol_algo == :lit # naive linear time iteration (deprecated)
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
            @timev begin
            while abs(diff1)>1e-6 && i<1000
                Mat[:,1:n_par.nstates] = BB[:,1:n_par.nstates] .+ CC * F0
                F00 = Mat \ (-AA)
                #F00 .= (BB .+ CC * F0) \ (-AA)
                diff1 = maximum(abs.(F00[:] .- F0[:]))[1]
                F0 .= F00
                i += 1
            end
            end
            hx = F0[1:n_par.nstates,1:n_par.nstates]
            gx = F0[n_par.nstates+1:end,1:n_par.nstates]
            println(i)
            println("lit")
            nk = n_par.nstates
            alarm_sgu = false
        end
    elseif n_par.sol_algo == :litx # linear time iteration with speed-up following Reiter
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
            X21                                      = zeros(n_par.ncontrols, n_par.nstates)
            X21up                                    = zeros(n_par.ncontrols, n_par.nstates)
            X21[1:n_par.ncontrols, 1:n_par.nstates] .= n_par.State2Control_save

            diff1 = 1000.0
            i = 0
            Q1 = Z[1:n_par.nstates,:]*F2        # [Z11 Z12] * F2
            Q2 = Z[n_par.nstates+1:end,:]*F2    # [Z21 Z22] * F2
            Q3 = Z[1:n_par.nstates,:]*F3        # [Z11 Z12]*F3
            Q4 = Z[n_par.nstates+1:end,:]*F3    # [Z21 Z22]*F3
            
            H  = I - X21*((I+Q1*X21)\Q1)
            # within loop only update X21
            while diff1>1e-6 && i<1000 
                i       += 1
                X21up   .= Q4 - Q2*H*X21*Q3     # [Z21 - Q2*H*X21*Z11  Z22 - Q2*H*X21*Z12]*F3
                diff1    = maximum(abs.(X21up .- X21))[1]
                X21     .= X21up
                H       .= I - X21*((I+Q1*X21)\Q1)
            end
            
            X11          = Q3 - Q1*H*X21*Q3     # [Z11 - Q1 * H*X21*Z11  Z12 - Q1*H*X21*Z12]* F3
            hx           = X11                  # F0[1:n_par.nstates,1:n_par.nstates]
            gx           = X21                  # [n_par.nstates+1:end,1:n_par.nstates]
            nk           = n_par.nstates
            alarm_sgu    = false
        end
    elseif n_par.sol_algo == :schur # (complex) schur decomposition
        alarm_sgu = false
        Schur_decomp, slt, nk, λ = complex_schur(A, -B) # first output is generalized Schur factorization

        # Check for determinacy and existence of solution
        if n_par.nstates != nk
            if estim # return zeros if not unique and determinate
                hx = Array{Float64}(undef, n_par.nstates, n_par.nstates)
                gx = Array{Float64}(undef,  n_par.ncontrols, n_par.nstates)
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
    return gx, hx, alarm_sgu, nk
end