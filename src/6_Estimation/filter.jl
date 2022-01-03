@doc raw"""
    likeli(par,Data,Data_miss,H_sel,XSSaggr,A,B,indexes,indexes_aggr,m_par,n_par,e_set,Copula,distrSS,compressionIndexes,priors,meas_error,meas_error_std;smoother=False)

Compute the likelihood of `Data`, given model-parameters `par` and prior `priors` (maximize to find MLE of `par`).

Solve model with [`SGU_estim()`](@ref), compute likelihood with [`kalman_filter()`](@ref) or with [`kalman_filter_smoother()`](@ref) (if `smoother==True`).

# Returns
*if `smoother==False`:*
- `log_like`,`prior_like`,`post_like`,`alarm`: log-likelihoods (`post` is the sum of `prior` and computed likelihood); `alarm` indicates error when solving model with [`SGU_estim`](@ref), sets log-likelihood to `-9.e15`
*if `smoother==True`:*
- `smoother_output`: returns from [`kalman_filter_smoother()`](@ref)
"""
function likeli(par, Data, Data_miss, H_sel, XSSaggr, A, B, indexes,indexes_aggr, m_par, n_par, e_set,distrSS, compressionIndexes, priors, meas_error, meas_error_std; smoother=false)

    # check priors, abort if they are violated
    prior_like::eltype(par), alarm_prior::Bool = prioreval(Tuple(par), Tuple(priors))
    alarm = false
    if alarm_prior
        log_like = -9.e15
        alarm = true
        State2Control = zeros(n_par.ncontrols, n_par.nstates)
        if e_set.debug_print
            println("PRIOR")
        end
    else
        if e_set.me_treatment != :fixed
            m_start = length(par) - length(meas_error) # find out where in par structural pars end
        else
            m_start = length(par)
        end

        # replace estimated values in m_par by last candidate
        m_par = Flatten.reconstruct(m_par, par[1:m_start])

        # covariance of structural shocks
        SCov = zeros(eltype(par), n_par.nstates, n_par.nstates)
        for i in e_set.shock_names
        	SCov[getfield(indexes, i), getfield(indexes, i)] = (getfield(m_par,Symbol("σ_", i))).^2
        end

        # covariance of measurement errors, assumption: ME ordered after everything else
        m = size(H_sel)[1]
        MCov = Diagonal(zeros(eltype(par), m)) # no correlated ME allowed for now
        if !isempty(meas_error)
            m_iter = 1
            if e_set.me_treatment != :fixed
                for (k, v) in meas_error # read out position of measurement errors
                    MCov[v, v] = par[m_start + m_iter].^2
                    m_iter += 1
                end
            else
                for (k, v) in meas_error # read out position of measurement errors
                    MCov[v, v] = meas_error_std[m_iter].^2
                    m_iter += 1
                end
            end
        end

        # solve model using candidate parameters
        # BLAS.set_num_threads(1)
        State2Control::Array{eltype(par),2}, LOM::Array{eltype(par),2}, alarm_sgu::Bool = SGU_estim(XSSaggr, A,B, m_par, n_par, indexes,
                                                                                            indexes_aggr, distrSS; estim = true)

        # BLAS.set_num_threads(Threads.nthreads())
        if alarm_sgu # abort if model doesn't solve
            log_like = -9.e15
            alarm = true
            if e_set.debug_print
                println("SGU")
            end
        else
            MX = [I; State2Control]
            H = H_sel * MX
            if smoother == false
                log_like = kalman_filter(H, LOM, Data, Data_miss, SCov, MCov, e_set)
                # log_like = kalman_filter_herbst(Data, LOM, SCov, H, MCov, 0, e_set)
            else
                smoother_output = kalman_filter_smoother(H, LOM, Data, .!Data_miss, SCov, MCov, e_set)
                log_like = smoother_output[1]
            end
        end
    end
    post_like = log_like + prior_like

    if smoother == false
        return log_like, prior_like, post_like, alarm, State2Control
    else
        return smoother_output
    end

end

@doc raw"""
    kalman_filter(H,LOM,Data,D_miss,SCov,MCov,e_set)

Compute likelihood of `Data`, applying the Kalman filter to the state-space represenation (`H`,`LOM`)
of the model.

# Arguments
- `H::Array{Float64,2}`: observation equation
- `LOM::Array{Float64,2}`: law of motion for states
- `Data::Array{Union{Missing,Float64},2}`,`D_miss::BitArray{2}`: data (time ``\times`` variable); marker for missing data
- `SCov::Array{Float64,2}`: covariance of structural shocks
- `MCov::Diagonal{Float64,Array{Float64,1}}`: covariance of measurement error

# Returns
- log-likelihood
"""
function kalman_filter(H::Array, LOM::Array, Data::Array{Union{Missing,Float64},2},
    D_miss::BitArray{2}, SCov::Array, MCov::Diagonal, e_set::EstimationSettings)

    # treat non-well-behaved covariance matrix
    SIG = lyapd(LOM, SCov)
    SIG = nearest_spd(SIG)

    t = size(Data)
    n = size(LOM)
    xhat = zeros(Float64, n[1])
    log_lik = 0.0

    H_slice = copy(H)
    @views @inbounds for s = 1:t[1]
        miss_temp = findall(D_miss[s, :])
        Data_slice = Data[s, :]
        Data_slice[miss_temp] .= 0.0
        copyto!(H_slice, H)
        H_slice[miss_temp, :] .= 0.0

        # compute likelihood contribution
        resi = Data_slice .- H_slice * xhat
        SH = SIG * H_slice'
        Ω = H_slice * SH + MCov

        for i in miss_temp
            Ω[i, :] .= 0.0
            Ω[:, i] .= 0.0
            Ω[i, i] = 1.0
        end
        OmegaInv = Ω \ I

        logdet_Ω, sign_logdet = logabsdet(Ω)
        if sign_logdet < 0
            log_lik += -10.e8
            if e_set.debug_print
                println("KF")
            end
            return log_lik
        else
            log_lik += -logdet_Ω - resi' * OmegaInv * resi - (t[2] - length(miss_temp)) * log(2.0 * π)
        end

        # update
        K = LOM * SH * OmegaInv # Gain
        xhat = LOM * xhat .+ K * resi
        Z = LOM .- K * H_slice
        SIG .= Z * (SIG * Z') .+ K * (MCov * K') .+ SCov
        SIG .= 0.5 .* (SIG .+ SIG')
    end
    loglik = 0.5 * log_lik

    return loglik
end

# kalman filter variant without missing data
function kalman_filter(H::Array, LOM::Array, Data::Array{Float64,2},
    D_nomiss::BitArray{2}, SCov::Array, MCov::Diagonal, e_set::EstimationSettings)

    # treat non-well-behaved covariance matrix
    SIG = lyapd(LOM, SCov)
    SIG = nearest_spd(SIG)

    t = size(Data)
    n = size(LOM)
    xhat = zeros(Float64, n[1])
    log_lik = 0.0

    @views @inbounds for s = 1:t[1]
        # compute likelihood contribution
        resi = Data[s, :] .- H * xhat
        SH   = SIG * H'
        Ω    = H * SH + MCov
        OmegaInv = I / Ω

        logdet_Ω, sign_logdet = logabsdet(Ω)
        if sign_logdet < 0
            log_lik += -10.e8
            if e_set.debug_print
                println("KF")
            end
            return log_lik
        else
            log_lik += -logdet_Ω - resi' * OmegaInv * resi
        end

        # update
        K = LOM * SH * OmegaInv # Gain
        xhat = LOM * xhat + K * resi
        Z = LOM - K * H
        SIG = Z * (SIG * Z') + K * (MCov * K') + SCov
        SIG = 0.5 * (SIG + SIG')
    end
    loglik = 0.5 * log_lik - 0.5 * t[1] * t[2] * log(2.0 * π)

    return loglik
end

@doc raw"""
    kalman_filter_smoother(H,LOM,Data,D_nomiss,SCov,MCov,e_set)

Compute likelihood and estimate of underlying states given the full observed `Data` by
applying the Kalman smoother to the state-space represenation (`H`,`LOM`) of the model.

# Arguments
- `H::Array{Float64,2}`: observation equation
- `LOM::Array{Float64,2}`: law of motion for states
- `Data::Array{Union{Missing,Float64},2}`,`D_nomiss::BitArray{2}`: data (time ``\times`` variable); marker for existent data
- `SCov::Array{Float64,2}`: covariance of structural shocks
- `MCov::Diagonal{Float64,Array{Float64,1}}`: covariance of measurement error

# Returns
- `log_lik`: log-likelihood
- `xhat_tgt`,`xhat_tgT`: estimate of underlying states from forward iteration [`xhat_tgt`] and
    backward iteration [`xhat_tgT`]
- `Sigma_tgt`,`Sigma_tgT`: estimate of covariance matrices from forward iteration [`Sigma_tgt`]
    and backward iteration [`Sigma_tgT`]
- `s`,`m`: ?
"""
function kalman_filter_smoother(H::Array{Float64,2}, LOM::Array{Float64,2}, Data, D_nomiss::BitArray{2},
                                SCov::Array{Float64,2}, MCov::Diagonal{Float64, Array{Float64,1}}, e_set)

    T, n_obs_vars = size(Data)
    n_states = size(LOM)[1]

    Sigma_tgtm1 = zeros(Float64, n_states, n_states, T+1)
    SIG::Array{Float64,2} = lyapd(LOM, SCov)
    Sigma_tgtm1[:, :, 1] = nearest_spd(SIG)

    Sigma_tgt = zeros(Float64, n_states, n_states, T)
    K = zeros(Float64, n_states, n_obs_vars, T)
    L = zeros(Float64, n_states, n_states, T)
    xhat_tgtm1 = zeros(Float64, n_states, T+1)
    xhat_tgt = zeros(Float64, n_states, T)
    resi = zeros(Float64, n_obs_vars, T)
    OmegaInv = zeros(Float64, n_obs_vars, n_obs_vars, T)
    log_lik = 0.0

    for t = 1:T
       # compute likelihood contribution
       resi[D_nomiss[t, :], t] = Data[t, D_nomiss[t, :]] .- H[D_nomiss[t, :], :] * xhat_tgtm1[:, t]
       SH   = Sigma_tgtm1[:, :, t] * transpose(H[D_nomiss[t, :], :])
       Ω    = H[D_nomiss[t, :], :] * SH + MCov[D_nomiss[t, :], D_nomiss[t, :]]
       logdet_Ω, sign_logdet = logabsdet(Ω)
       if sign_logdet < 0
            log_lik+=-10.e8
            if e_set.debug_print
                println("KF")
            end
            return log_lik
        else
            OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] = I / Ω
            log_lik += - 0.5 * n_obs_vars * log(2.0 * π) -
                        0.5 * logdet_Ω - 0.5 * transpose(resi[D_nomiss[t, :], t]) *
                        OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] * resi[D_nomiss[t, :], t]
        end

       # update
       K[:, D_nomiss[t, :], t]  = LOM * SH * OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] # Gain
       L[:, :, t]  = LOM - K[:, D_nomiss[t, :], t] * H[D_nomiss[t, :], :]
       xhat_tgt[:, t] = xhat_tgtm1[:, t] + (SH * OmegaInv[D_nomiss[t, :], D_nomiss[t, :],  t]) * resi[D_nomiss[t, :], t]
       Sigma_tgt[:, :, t] = Sigma_tgtm1[:, :, t] - SH * (OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] * transpose(SH))
       xhat_tgtm1[:, t + 1] = LOM * xhat_tgtm1[:, t] + K[:, D_nomiss[t, :], t] * resi[D_nomiss[t, :], t]
       Sigma_tgtm1_temp =  L[:, :, t] * (Sigma_tgtm1[:, :, t] * L[:, :, t]') + K[:, D_nomiss[t, :], t] *
                            (MCov[D_nomiss[t, :], D_nomiss[t, :]] * K[:, D_nomiss[t, :], t]') + SCov
       Sigma_tgtm1[:, :, t + 1] = 0.5 * (Sigma_tgtm1_temp + transpose(Sigma_tgtm1_temp))

    end

    xhat_tgT = zeros(Float64, n_states, T)
    xhat_tgT[:, T] = xhat_tgt[:, T]
    Sigma_tgT = zeros(Float64, n_states, n_states, T)
    Sigma_tgT[:, :, T] = Sigma_tgt[:, :, T]
    r = zeros(Float64, n_states, T + 1)
    N = zeros(Float64, n_states, n_states, T+1)
    s = zeros(Float64, n_states, T)
    m = zeros(Float64, n_obs_vars, T)
    for t = T:-1:1
       r[:, t] = H[D_nomiss[t, :], :]' * OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] * resi[D_nomiss[t, :], t] + L[:, :, t]' * r[:, t + 1]
       xhat_tgT[:, t] = xhat_tgtm1[:, t] + Sigma_tgtm1[:, :, t] * r[:, t]

       N[:, :, t] = H[D_nomiss[t, :], :]' * OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] * H[D_nomiss[t, :], :] + L[:, :, t]' * N[:, :, t + 1] * L[:, :, t]
       Sigma_tgT[:, :, t] = Sigma_tgtm1[:, :, t] - Sigma_tgtm1[:, :, t] * N[:, :, t] * Sigma_tgtm1[:, :, t]

       s[:, t] = SCov * r[:, t+1]
       m[D_nomiss[t, :], t] = MCov[D_nomiss[t, :], D_nomiss[t, :]] * (OmegaInv[D_nomiss[t, :], D_nomiss[t, :], t] * resi[D_nomiss[t, :], t] - K[:, D_nomiss[t, :], t]' * r[:, t+1])
    end

    return log_lik, xhat_tgt, xhat_tgT, Sigma_tgt, Sigma_tgT, s, m
end

function kalman_filter_herbst(Data, LOM, SCov, H, MCov, t0, e_set)
    tol = 1e-7
    converged = false

    nobs, ny = size(Data)
    ns = size(LOM, 1)
    xhat = zeros(ns, 1)

    P::Array{Float64,2} = lyapd(LOM, SCov)

    F = H * (P * H') + MCov
    F = 0.5 * (F + F')
    iF = inv(F)

    K = LOM * (P * H')
    W = copy(K)
    M = -iF
    Kg = K * iF

    log_lik = 0.0
    for i = 1:nobs

        # calculate the forecast errors
        ν = Data[i, :] .- H * xhat

        logdet_F, sign_logdet = logabsdet(F)
        if sign_logdet < 0
             log_lik += -10.e8
             if e_set.debug_print
                 println("KF")
             end
             return log_lik
         else
             iFν = F \ ν
             if i > t0
                 log_lik += - logdet_F - (ν' * iFν)[1]
             end
         end

        # updating
        xhat = LOM * xhat + Kg * ν

        if !converged
            # these are intermediate calculations we can re-use
            HWM = H * W * M
            HWMW = HWM * W'

            M = M + (HWM' * iF) * HWM  # M_[t+1]
            M = 0.5 * (M + M')
            F = F + HWMW * H'         # F_[t+1]
            F = 0.5 * (F + F')
            iF = inv(F)
            K = K + LOM * HWMW'       # K_[t+1]
            Kg_old = copy(Kg)
            Kg = K * iF
            W = (LOM - Kg * H) * W    # W_[t+1]

            if maximum(abs.(Kg .- Kg_old)) < tol
                converged = true
            end
        end
    end
    log_lik = 0.5 * log_lik - 0.5 * nobs * ny * log(2 * pi)
    return log_lik
end

@doc raw"""
    nearest_spd(A)

Return the nearest (in Frobenius norm) Symmetric Positive Definite matrix to `A`.

Based on [answer in MATLAB Central forum](https://de.mathworks.com/matlabcentral/answers/320134-make-sample-covariance-correlation-matrix-positive-definite).
From Higham: "The nearest symmetric positive semidefinite matrix in the
Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
where H is the symmetric polar factor of B=(A + A')/2."

# Arguments
`A`: square matrix

# Returns
`Ahat`: nearest SPD matrix to `A`
"""
function nearest_spd(A)

    # symmetrize A into B
    B = 0.5 .* (A .+ A')
    FU, FS, FVt = LinearAlgebra.LAPACK.gesvd!('N', 'S', copy(B))
    H = FVt' * Diagonal(FS) * FVt

    # get Ahat in the above formula
    Ahat = 0.5 .* (B .+ H)

    # ensure symmetry
    Ahat .= 0.5 .* (Ahat .+ Ahat')
    
    # test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
    p = false
    k = 0
    count = 1
    while p == false && count<100
        R = cholesky(Ahat; check=false)
        k += 1
        count = count + 1
        if ~issuccess(R)
            # Ahat failed the chol test. It must have been just a hair off,
            # due to floating point trash, so it is simplest now just to
            # tweak by adding a tiny multiple of an identity matrix.
            mineig = eigmin(Ahat)
            Ahat .+= (-mineig .* k.^2 .+ eps(mineig)) .* Diagonal(ones(size(Ahat,1)))
        else
            p = true
        end
    end

    return Ahat
end