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
function likeli(par, Data, Data_miss, H_sel, XSSaggr, A, B, indexes,indexes_aggr, m_par, n_par, e_set, Copula,distrSS, compressionIndexes, priors, meas_error, meas_error_std; smoother=false)
    # display("-----------------------")
    # display("Parameters")
    # display(par')
    prior_like::Float64, alarm_prior::Bool = prioreval(Tuple(par), Tuple(priors))
    # write par to model parameters
    alarm = false

    if alarm_prior
        log_like = -9.e15
        alarm = true

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

        #display(m_par)
        # Covariance of structural shocks
        SCov = zeros(n_par.nstates, n_par.nstates)
        for i in e_set.shock_names
        	SCov[getfield(indexes, i), getfield(indexes, i)] = (getfield(m_par,Symbol("σ_", i))).^2
        end

        # Covariance of measurement errors, assumption: ME ordered after everything else
        m = size(H_sel)[1]
        MCov = Diagonal(zeros(m)) # no correlated ME allowed for now
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
        # solve this
        BLAS.set_num_threads(1)
        State2Control::Array{Float64,2}, LOM::Array{Float64,2}, alarm_sgu::Bool = SGU_estim(XSSaggr, A,B, m_par, n_par, indexes, indexes_aggr, distrSS; estim = true)

        BLAS.set_num_threads(Threads.nthreads())

        if alarm_sgu
            log_like = -9.e15
            alarm = true

            if e_set.debug_print
                println("SGU")
            end
        else
            MX = [I; State2Control] # States t; Controls t
            if e_set.fd_flag
                LOM_D  = copy([LOM zeros(size(LOM)); I zeros(size(LOM))])
                SCov_D = [SCov zeros(size(SCov)); zeros(size(SCov)) zeros(size(SCov))]
                MX_D   = copy([MX -MX; MX zeros(size(MX))])
                H = H_sel * MX_D
                if smoother == false
                    log_like = kalman_filter(H, LOM_D, Data, Data_miss, SCov_D, MCov, e_set)
                    # log_like = kalman_filter_herbst(Data, LOM_D, SCov_D, H, MCov, 0, e_set)
                else
                    smoother_output = kalman_filter_smoother(H, LOM_D, Data, .!Data_miss, SCov_D, MCov, e_set)
                    log_like = smoother_output[1]
                end
            else
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
    end
    post_like = log_like + prior_like

    if smoother == false
        return log_like, prior_like, post_like, alarm
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
function kalman_filter(H::Array{Float64,2}, LOM::Array{Float64,2}, Data::Array{Union{Missing, Float64},2},
    D_miss::BitArray{2}, SCov::Array{Float64,2}, MCov::Diagonal{Float64, Array{Float64,1}}, e_set::EstimationSettings)

    SIG = lyapunov_symm_own(LOM, SCov, 1e-15)
    SIG = nearest_spd(SIG)
    # display(isposdef(SIG))

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
       SH   = SIG * H_slice'
       Ω    = H_slice * SH + MCov

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
            log_lik += - logdet_Ω - resi' * OmegaInv * resi - (t[2] - length(miss_temp)) * log(2.0 * π)
        end
       # update
       K  = LOM * SH * OmegaInv# Gain
       xhat = LOM * xhat .+ K * resi
       Z = LOM .- K * H_slice
       SIG .= Z * (SIG * Z') .+ K * (MCov * K')  .+ SCov
       SIG .= 0.5 .* (SIG .+ SIG')
    end
    loglik = 0.5 * log_lik #- 0.5 * t[1] * t[2] * log(2.0 * π)

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
function kalman_filter_smoother(H::Array{Float64,2}, LOM::Array{Float64,2},
    Data::Array{Union{Missing, Float64},2}, D_nomiss::BitArray{2}, SCov::Array{Float64,2},
    MCov::Diagonal{Float64, Array{Float64,1}}, e_set)

    T, n_obs_vars = size(Data)
    n_states = size(LOM)[1]

    Sigma_tgtm1 = zeros(Float64, n_states, n_states, T+1)
    SIG::Array{Float64,2} = lyapunov_symm_own(LOM,SCov,1e-15)
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
       # if t>100
       #     resi[end-4,t]=0
       # end
       # println(t)
       # println(resi[1, t])
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
       # Sigma_tgtm1[:, :, t + 1] = LOM * (Sigma_tgt[:, :, t] * LOM') + SCov # DK 4.23
       # Sigma_tgtm1_temp = LOM * (Sigma_tgt[:, :, t] * transpose(LOM)) + SCov # DK 4.23
       Sigma_tgtm1_temp =  L[:, :, t] * (Sigma_tgtm1[:, :, t] * L[:, :, t]') + K[:, D_nomiss[t, :], t] * (MCov[D_nomiss[t, :], D_nomiss[t, :]] * K[:, D_nomiss[t, :], t]') + SCov
       # Sigma_tgtm1[:, :, t + 1] =  L[:, :, t] * (Sigma_tgtm1[:, :, t] * L[:, :, t]') + K[:, :, t] * (MCov * K[:, :, t]') + SCov
       # Sigma_tgtm1[:, :, t + 1] = tril(Sigma_tgtm1_temp) + transpose(tril(Sigma_tgtm1_temp, -1))
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

    P::Array{Float64,2} = lyapunov_symm_own(LOM, SCov, 1e-15)

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

function lyapunov_symm_own(a::Array{Float64,2}, b::Array{Float64,2}, lyapunov_complex_threshold::Float64)

    T, Z = schur(a)
    B = (Z' * b) * Z
    n = size(T, 1)
    x = zeros(n, n)
    i = copy(n)

    @views @inbounds while i >= 2
        if abs(T[i, i-1]) < lyapunov_complex_threshold
            if i .== n
                c = zeros(n, 1)
            else
                c = T[1:i, :] * (x[:, i+1:end] * T[i,i+1:end]) +
                    T[i, i] .* (T[1:i, i+1:end] * x[i+1:end,i])
            end
            q           = I - T[1:i, 1:i] .* T[i, i]
            x[1:i, i]   = q \ (B[1:i, i] .+ c)
            x[i, 1:i-1] = x[1:i-1, i]'
            i -= 1
        else
            if i .== n
                c   = zeros(n, 1)
                c1  = zeros(n, 1)
            else
                c   = T[1:i, :] * (x[:, i+1:end] * T[i, i+1:end]) +
                        T[i, i] .* (T[1:i, i+1:end] * x[i+1:end, i]) +
                        T[i, i-1] .* (T[1:i, i+1:end] * x[i+1:end, i-1])
                c1  = T[1:i, :] * (x[:, i+1:end] * T[i-1, i+1:end]) +
                        T[i-1, i-1] .* (T[1:i, i+1:end] * x[i+1:end, i-1]) +
                        T[i-1, i] .* (T[1:i, i+1:end] * x[i+1:end, i])
            end
            q = Matrix{Float64}(undef, 2*i, 2*i)
            for i2 = 1:i, i1 = 1:i
                if i1 == i2
                    q[i1, i2]     = 1.0 .- T[i, i] .* T[i1, i2]
                    q[i1, i+i2]   = - T[i, i-1] .* T[i1, i2]
                    q[i+i1, i2]   = - T[i-1, i] .* T[i1, i2]
                    q[i+i1, i+i2] = 1.0 .- T[i-1, i-1] .* T[i1, i2]
                else
                    q[i1, i2]     = - T[i, i] .* T[i1, i2]
                    q[i1, i+i2]   = - T[i, i-1] .* T[i1, i2]
                    q[i+i1, i2]   = - T[i-1, i] .* T[i1, i2]
                    q[i+i1, i+i2] = - T[i-1, i-1] .* T[i1, i2]
                end
            end
            z =  q \ [ B[1:i, i] .+ c; B[1:i, i-1] .+ c1 ]
            x[1:i, i]       .= z[1:i]
            x[1:i, i-1]     .= z[i+1:end]
            x[i, 1:i-1]     .= x[1:i-1, i]
            x[i-1, 1:i-2]   .= x[1:i-2, i-1]
            i -= 2
        end
    end
    if i .== 1
        c = dot(T[1, :], x[:, 2:end] * T[1, 2:end]) .+ T[1, 1] .* dot(T[1, 2:end], x[2:end, 1])
        x[1, 1] = (B[1, 1] .+ c) ./ (1 .- T[1, 1] .* T[1, 1])
    end
    x = Z * x * Z'

    return x
end
