@doc raw"""
    SolveDiffEq(A, B, n_par; estim)

Calculate the solution to the linearized difference equations defined
B x_t = A x_{t+1}.

# Arguments
- `A`,`B`: matrices with first derivatives 
- `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm, options are: 
    * `litx`:  Linear time iteration (improved implementation following Reiter fast if initial guess is good)
    * `schur`: Klein's algorithm (fastest if no good initial guess is available)

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_sgu`,`nk`: `alarm_sgu=true` when solving algorithm fails, `nk` number of
    predetermined variables
"""
function SolveDiffEq(A::Array,B::Array, n_par::NumericalParameters, estim=false)
    if n_par.sol_algo == :litx # linear time iteration with speed-up following Reiter
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

            F1 = A[:, 1:n_par.nstates]
            F2 = A[:, n_par.nstates+1:end]
            F3 = -B[:, 1:n_par.nstates]
            F4 = B[:, n_par.nstates+1:end]
    
            Z = hcat(F1, F4) \ I
            Z1 = Z[1:n_par.nstates, :]
            Z2 = Z[n_par.nstates+1:end, :]
            Q1 = Z1 * F2    # [Z11 Z12] * F2
            Q2 = Z2 * F2    # [Z21 Z22] * F2
            Q3 = Z1 * F3    # [Z11 Z12] * F3
            Q4 = Z2 * F3    # [Z21 Z22] * F3

    
            diff1 = 1000.0
            i = 0

            # write model as Z = (X2 - Q4) = - (Q2*Q4 + Q2*(X2-Q4))*X11
            # use QR decomposition on [Q2Q4 Q2] to find which X2-Q4 are
            # linearly dependent for sure and by premultiplying 
            # Q' can be eliminated from the system:
            # Q'Z = -(R1 + R2*Q*Z)*X11
            # This not only makes the system smaller, but also increases precision.
            # Then make use of the possibility to fuse left hand multiplications.  
            F       = qr([Q2*Q4 Q2], ColumnNorm())
            red     = count(abs.(diag(F.R)).<1.0e-9)
            Q       = Matrix(F.Q)
            R       = F.R*F.P'
            R1      = R[1:end-red,1:n_par.nstates] # Split R matrix 
            R2      = R[1:end-red,n_par.nstates+1:end] * Q[:,1:end-red]
            Q1hat   = Q1*Q[:,1:end-red] # Matrix to build back S2S
            Q14     = I+Q1*Q4
            Z       = -R1*(Q14 \ Q3) # initial guess
            X11     = (Q1hat * Z + Q14) \ Q3
            
            R21     = R2*R1
            R221    = (R2^2)*R1
            R2221   = (R2^3)*R1
            R22221  = (R2^4)*R1
            R22222  = R2^5
        
            # stacked iteraded (-1)^tR1*R2^t 
            T       =-[R1 -R21 R221 -R2221 R22221]
            F2      = qr(T, ColumnNorm()) # T has massive reduced rank (due to economic structure)
            red2    = count(abs.(diag(F2.R)).<1.0e-9) # find size of nullspace
            T2      = (F2.R*F2.P')[1:end-red2,:] #keep only Rows in traingular form that correspond to range
            QQ      = Matrix(F2.Q)[:,1:end-red2] # keep only columns of F2.Q that correspond to range, (others are 0 in T2)
            # Stacked iterated LOM for 4 periods
            X1c     = X11^5
            XX      =  [X11; X11^2;X11^3;X11^4;  X1c]
            while diff1 > 1e-11 && i < 500
                i += 1
                Zc = QQ*(T2*XX)
                Z  = Zc - R22222*Z*X1c#- (R1 - (R21 - (R221 -(R2221 + R2222 * Z)*X11) *X11) * X11) * X11
            
                for t = 1:min(5, ceil(1  - 0.5*log10(diff1))) # howard's improvement
                    Z = Zc - R22222*Z*X1c
                end
            
                up = (Q1hat * Z + Q14) \ Q3
                diff1 = maximum(abs.(up .- X11))[1]./maximum(abs.(up))[1]
                X11 .=  up
                X1c .= X11^5
                XX  .=  [X11; X11^2;X11^3;X11^4;  X1c]
                
                # Making use of 
                # H*X21 = (I - X21*(inv(I+Q1*X21)*Q1)*X21 
                #       = X21 - X21*inv(I+Q1*X21)*Q1*X21  - X21*inv(I+Q1*X21)) + X21*inv(I+Q1*X21)
                #             = X21 - X21*inv(I+Q1*X21)(I+Q1*X21) + X21*(inv(I+Q1*X21)))
                #             = X21 - X21 - X21*inv(I+Q1*X21) = X21*inv(I+Q1*X21)
                # X11= Q3 - Q1*X21*inv(I+Q1*X21) *Q3 = Q3 - (I-inv(I+Q1*X21)) *Q3 = inv(I+Q1*X21) *Q3
            end
            hx = X11
            gx = Q[:,1:end-red]*Z +Q4
            nk = n_par.nstates
            alarm_sgu = false
            if any(isnan.(X11))||any(isinf.(X11))
                nk = n_par.nstates-1
                alarm_sgu = true
                println("divergence of X11 to infty")
            elseif maximum(abs.(eigvals(X11)))>1
                nk = n_par.nstates-1
                alarm_sgu = true
                println("No stable solution")
                println(maximum(abs.(eigvals(X11))))
            elseif i==350
                alarm_sgu = true
                println("LITX not converged")
            end
            # println(i)
        end
    elseif n_par.sol_algo == :schur # (complex) schur decomposition
        alarm_sgu = false
        Schur_decomp, slt, nk, λ = complex_schur(A, -B) # first output is generalized Schur factorization
    
        # Check for determinacy and existence of solution
        if n_par.nstates != nk
            if estim # return zeros if not unique and determinate
                hx = zeros(eltype(A), n_par.nstates, n_par.nstates)
                gx = zeros(eltype(A), n_par.ncontrols, n_par.nstates)
                alarm_sgu = true
                return gx, hx, alarm_sgu, nk, A, B
            else # debug mode/ allow IRFs to be produced for roughly determinate system
                ind = sortperm(abs.(λ); rev = true)
                slt = zeros(Bool, size(slt))
                slt[ind[1:n_par.nstates]] .= true
                alarm_sgu = true
                @warn "critical eigenvalue moved to:"
                print(λ[ind[n_par.nstates-5:n_par.nstates+5]])
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
            hx = zeros(eltype(A), n_par.nstates, n_par.nstates)
            gx = zeros(eltype(A), n_par.ncontrols, n_par.nstates)
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