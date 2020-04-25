@doc raw"""
    SGU(XSS,A,B,m_par,n_par,indexes,Copula,compressionIndexes,distrSS;estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys()`](@ref), using Schmitt-Grohé & Uribe (JEDC 2004) style linearization
(apply the implicit function theorem to obtain linear observation and
state transition equations).

The Jacobian is calculated using dual numbers (implemented by package `ForwardDiff`).
Use macro `@make_deriv` to compute partials simultaneously, with chunk size
given by `global` `n_FD`. Make use of model knowledge to set some entries manually.

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

    ############################################################################
    # Check whether Steady state solves the difference equation
    ############################################################################
    length_X0 = indexes.profits # Convention is that profits is the last control
    X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,tuple(zeros(n_FD)...))
    @make_deriv n_FD
    F  = Fsys(X0,X0,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC, Copula)
    # if maximum(abs.(F))/10>n_par.ϵ
    #     @warn  "F=0 is not at required precision"
    # end
    BLAS.set_num_threads(1)


    ############################################################################
    # Calculate Jacobians of the Difference equation F
    ############################################################################
    # B = zeros(length_X0, length_X0)
    # A = zeros(length_X0, length_X0)
    #-------------- loop over states T -----------------------------------------
    state_loop!(B, Deriv, indexes,n_FD)

    #------------- loop over controls T ----------------------------------------

    # Make use of the fact that Vk/Vm has no influence on any variable in
    # the system, thus derivative is 1
    for i in indexes.Vm
        B[i, i] = 1.0
    end

    for i in indexes.Vk
        B[i, i] = 1.0
    end

    control_loop!(B, Deriv, indexes, length_X0,n_FD)

    # ------------ loop over states and controls T+1 ---------------------------
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
    prime_loop!(A, DerivPrime, indexes, length_X0,n_FD)

    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    BLAS.set_num_threads(Threads.nthreads())
    # BLAS.set_num_threads(1)
    if n_par.sol_algo == :lit
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
            while abs(diff1)>1e-6 && i<1000
                Mat[:,1:n_par.nstates] = BB[:,1:n_par.nstates] .+ CC * F0
                F00 = Mat \ (-AA)
                #F00 .= (BB .+ CC * F0) \ (-AA)
                diff1 = maximum(abs.(F00[:] .- F0[:]))[1]
                F0 .= F00
                i += 1
            end
            hx = F0[1:n_par.nstates,1:n_par.nstates]
            gx = F0[n_par.nstates+1:end,1:n_par.nstates]

            nk = n_par.nstates
            alarm_sgu = false
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
        z21 = view(Schur_decomp.Z, (nk+1):length_X0, 1:nk)
        z11 = view(Schur_decomp.Z, 1:nk, 1:nk)
        s11 = view(Schur_decomp.S, 1:nk, 1:nk)
        t11 = view(Schur_decomp.T, 1:nk, 1:nk)

        if rank(z11) < nk
            @warn "invertibility condition violated"
            hx = Array{Float64}(undef, n_par.nstates, n_par.nstates)
            gx = Array{Float64}(undef, length(compressionIndexes) + 13, n_par.nstates)
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

function state_loop!(B, Deriv, indexes,n_FD)
    Threads.@threads for i = 1:n_FD:(indexes.Vm[1]-1) #in eachcol(B[:,1:(indexes.Vm[1]-1)])
        aux = Deriv.(i)
        for k = 1:min(n_FD,(indexes.Vm[1]-1)-i+1)
            for j = 1:size(B,1)
                B[j,i+k-1] = aux[j][k]
            end
        end
    end
end

function control_loop!(B, Deriv, indexes, length_X0,n_FD)
    Threads.@threads for i = (indexes.Vk[end]+1):n_FD:length_X0
        aux = Deriv.(i)
        for k = 1:min(n_FD,length_X0-i+1)
            for j = 1:size(B,1)
                B[j,i+k-1] = aux[j][k]
            end
        end
    end
end

# Calculating Jacobian to XPrime
function prime_loop!(A, DerivPrime, indexes, length_X0,n_FD)
    Threads.@threads for i = (1+indexes.distr_y[end]):n_FD:length_X0#in eachcol(A)
        aux = DerivPrime.(i)
        for k = 1:min(n_FD,length_X0-i+1)
            for j = 1:size(A,1)
                A[j,i+k-1] = aux[j][k]
            end
        end
    end
end
