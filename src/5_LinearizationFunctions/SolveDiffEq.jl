@doc raw"""
    SolveDiffEq(A, B, n_par; estim)

Calculate the solution to the linearized difference equations defined as
P'*B*P x_t = P'*A*P x_{t+1}, where `P` is the (ntotal x r) semi-unitary model reduction matrix
`n_par.PRightAll` of potentially reduced rank r.

# Arguments
- `A`,`B`: matrices with first derivatives 
- `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm, options are: 
    * `litx`:  Linear time iteration (implementation follows Reiter)
    * `schur`: Klein's algorithm (preferable if number of controls is small)

# Returns
- `gx`,`hx`: observation equations [`gx`] and state transition equations [`hx`]
- `alarm_sgu`,`nk`: `alarm_sgu=true` when solving algorithm fails, `nk` number of
    predetermined variables
"""
function SolveDiffEq(Ainput::Array,Binput::Array, n_par::NumericalParameters, estim=false)
    lit_fail = false
    A = n_par.PRightAll' * Ainput * n_par.PRightAll 
    B = n_par.PRightAll' * Binput * n_par.PRightAll
    if n_par.sol_algo == :lit # linear time iteration with speed-up following Reiter (less numerically stable than schur)
        @views begin
            # Formula X_{t+1} = [X11 X12;X21 X22] = inv(BB + [0 F2]*X_t)*[F3 0]
            # implies X21=0, X22 = 0
            # use of Woodburry-Morrison-Shermann
            # X.1_{t+1} =  (inv(BB)  - inv(BB)*[0 F2]* inv(I+[X.1_t 0]*inv(BB)*[0 F2])[X.1_t 0])*F3inv(BB) ; WMS formula
            # Z = inv(BB)
            # X.1_{t+1} =  Z*F3  - Z * [0 F2] * inv(I+[X.1_t 0]*Z*[0 F2])* [X.1_t 0] * Z * F3 ; multiply
            #           =  Z*F3  - Z * [0 F2] * inv(I+[0 X11_t *Z1. * F2;
            #                                          0 X21_t* Z1. * F2] ) * [X.1_t 0] * Z * F3
            #           =  Z*F3  - Z * [0 F2] * [I   -(X11_t * Z1. * F2) *inv(I+ X21_t* Z1. *F2); 
            #                                    0    inv(I+ X21_t* Z1. *F2)]) * [X.1_t 0] * Z * F3 ; block upper trianglular inverse
            #           =  Z*F3  - [0  , Z * F2 * inv(I+ X21_t* Z1. *F2)] * X.1_t * Z1. *F3 ; multiply
            #           =  Z*F3  - Z*F2* inv(I+ X21_t* Z1. *F2)* X21_t * Z1. *F3 ; multiply out
            #           =  [(Z1. * F3); (Z2. * F3) ] -  [(Z1. * F2); (Z2. * F2)]* (I - X21_t * inv(I +  (Z1. * F2) * X21_t)* (Z1. * F2)) * X21_t * (Z1. *F3); WMS formula
            #           =  [Q3; Q4 ] -  [Q1; Q2]* (I - X21_t * inv(I +  Q1 * X21_t)* Q1) * X21_t *Q3; WMS formula
            #
            # use inv(I+A)*A = I - inv(I+A)
            #
            # X11_{t+1} =  Q3 -  Q1 *(I - X21_t*inv(I +  Q1*X21_t) * Q1)*X21_t * Q3; 
            #           =  Q3 -  Q1 *(X21_t - X21_t*inv(I +  Q1*X21_t) * Q1*X21_t) * Q3;
            #           =  Q3 -  Q1 *(X21_t* inv(I+Q1*X21_t)) *Q3 
            #           =  Q3 -  (I - inv(I+Q1*X21_t)) *Q3 
            #           =  inv(I+Q1*X21_t) *Q3 
            # 
            # X21_{t+1} =  Q4 -  Q2 *(X21_t - X21_t*inv(I +  Q1*X21_t) * Q1*X21_t) * Q3; WMS formula
            #           =  Q4 -  Q2 *X21_t*inv(I +  Q1*X21_t) * Q3;  
      

            F1 = A[:, 1:n_par.nstates_r]
            F2 = A[:, n_par.nstates_r+1:end]
            F3 = -B[:, 1:n_par.nstates_r]
            F4 = B[:, n_par.nstates_r+1:end]
    
            Z = hcat(F1, F4) \ I
            Z1 = Z[1:n_par.nstates_r, :]
            Z2 = Z[n_par.nstates_r+1:end, :]
            Q1 = Z1 * F2    # [Z11 Z12] * F2
            Q2 = Z2 * F2    # [Z21 Z22] * F2
            Q3 = Z1 * F3    # [Z11 Z12] * F3
            Q4 = Z2 * F3    # [Z21 Z22] * F3

    
            diff1 = 1000.0
            i = 0

            X21     = Q4 # initial guess
            X11     = (I + Q1 *X21) \ Q3

            while diff1 > 1e-11 && i < 1000
                i += 1
                X21  = Q4 - Q2*X21 * X11           
                up = (I + Q1 *X21) \ Q3
                diff1 = maximum(abs.(up .- X11))[1]./maximum(abs.(up))[1]
                X11 .=  up                
            end
            hx = X11
            gx = X21
            nk = copy(n_par.nstates_r)
        end
        alarm_sgu = false
        if any(isnan.(X11))||any(isinf.(X11))
            nk = n_par.nstates_r-1
            lit_fail = true
            println("divergence of X11 to infty")
        elseif i==1000
            lit_fail = true
            println("LITX not converged -> trying with Schur")
        elseif maximum(abs.(eigvals(X11)))>1
            nk = n_par.nstates_r-1
            alarm_sgu = true
            println("No stable solution")
            println(maximum(abs.(eigvals(X11))))
        end
    end
    if n_par.sol_algo == :litx # linear time iteration with speed-up following Reiter, using QR and Howard-style improvement for speed-up (less stable)
        @views begin
            # Same setup as :lit
            F1 = A[:, 1:n_par.nstates_r]
            F2 = A[:, n_par.nstates_r+1:end]
            F3 = -B[:, 1:n_par.nstates_r]
            F4 = B[:, n_par.nstates_r+1:end]
    
            Z = hcat(F1, F4) \ I
            Z1 = Z[1:n_par.nstates_r, :]
            Z2 = Z[n_par.nstates_r+1:end, :]
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
            R1      = R[1:end-red,1:n_par.nstates_r] # Split R matrix 
            R2      = R[1:end-red,n_par.nstates_r+1:end] * Q[:,1:end-red]
            Q1hat   = Q1*Q[:,1:end-red] # Matrix to build back S2S
            Q14     = I+Q1*Q4
            Z       = -R1*(Q14 \ Q3) # initial guess
            X11     = (Q1hat * Z + Q14) \ Q3
            # Perform 3 update steps at once:
            R21     = R2*R1
            R221    = (R2^2)*R1
            R222    = R2^3
        
            # stacked iteraded (-1)^tR1*R2^t 
            T       =-[R1 -R21 R221 ]
            F2      = qr(T, ColumnNorm()) # T has massive reduced rank (due to economic structure)
            red2    = count(abs.(diag(F2.R)).<1.0e-9) # find size of nullspace
            T2      = (F2.R*F2.P')[1:end-red2,:] #keep only Rows in traingular form that correspond to range
            QQ      = Matrix(F2.Q)[:,1:end-red2] # keep only columns of F2.Q that correspond to range, (others are 0 in T2)
            # Stacked iterated LOM for 3 periods
            X1c     = X11^3
            XX      = [X11; X11^2; X1c]
            while diff1 > 1e-11 && i < 500
                i += 1
                Zc = QQ*(T2*XX)
                Z  = Zc - R222*Z*X1c#- (R1 - (R21 - (R221  * Z)*X11) *X11)  
                for t = 1:ceil(-8  - log10(diff1)) # if close to convergence: more update steps
                    Z = Zc - R222*Z*X1c
                end
                up = (Q1hat * Z + Q14) \ Q3
                diff1 = maximum(abs.(up .- X11))[1]./maximum(abs.(up))[1]
                X11 .=  up
                X1c .= X11^3
                XX  .=  [X11; X11^2; X1c]
            end
            hx = X11
            gx = Q[:,1:end-red]*Z +Q4
            nk = copy(n_par.nstates_r)
        end
        alarm_sgu = false
        if any(isnan.(X11))||any(isinf.(X11))
            nk = n_par.nstates_r-1
            lit_fail = true
            println("divergence of X11 to infty") 
        elseif i==500
            lit_fail = true
            println("LITX not converged -> trying with Schur")
        elseif maximum(abs.(eigvals(X11)))>1
            nk = n_par.nstates_r-1
            alarm_sgu = true
            println("No stable solution")
            println(maximum(abs.(eigvals(X11))))
        end
    end
    if n_par.sol_algo == :schur || lit_fail  # (complex) schur decomposition
        alarm_sgu = false
        Schur_decomp, slt, nk, λ = complex_schur(A, -B) # first output is generalized Schur factorization
        # Check for determinacy and existence of solution
        if n_par.nstates_r != nk
            if estim # return zeros if not unique and determinate
                hx = zeros(eltype(A), n_par.nstates_r, n_par.nstates_r)
                gx = zeros(eltype(A), n_par.ncontrols_r, n_par.nstates_r)
                alarm_sgu = true
                return gx, hx, alarm_sgu, nk, A, B
            else # debug mode/ allow IRFs to be produced for roughly determinate system
                ind = sortperm(abs.(λ); rev = true)
                slt = zeros(Bool, size(slt))
                slt[ind[1:n_par.nstates_r]] .= true
                alarm_sgu = true
                @warn "critical eigenvalue moved to:"
                print(λ[ind[n_par.nstates_r-5:n_par.nstates_r+5]])
                print(λ[ind[1]])
                nk = n_par.nstates_r
            end
        end
        # in-place reordering of eigenvalues for decomposition
        ordschur!(Schur_decomp, slt)
    
        # view removes allocations
        z21 = view(Schur_decomp.Z, (nk+1):n_par.ntotal_r, 1:nk)
        z11 = view(Schur_decomp.Z, 1:nk, 1:nk)
        s11 = view(Schur_decomp.S, 1:nk, 1:nk)
        t11 = view(Schur_decomp.T, 1:nk, 1:nk)
    
        if rank(z11) < nk
            @warn "invertibility condition violated"
            hx = zeros(eltype(A), n_par.nstates_r, n_par.nstates_r)
            gx = zeros(eltype(A), n_par.ncontrols_r, n_par.nstates_r)
            alarm_sgu = true
            return gx, hx, alarm_sgu, nk, A, B
        end
        z11i = z11 \ I # I is the identity matrix -> doesn't allocate an array!
        gx = real(z21 * z11i)
        hx = real(z11 * (s11 \ t11) * z11i)
    end
    if n_par.sol_algo != :schur && n_par.sol_algo != :litx && n_par.sol_algo != :litx_s
        error("Solution algorithm not defined!")
    end

    return gx, hx, alarm_sgu, nk
end