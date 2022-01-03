@doc raw"""
    rwmh(xhat, Σ, n_par, Data, D_miss, H_sel, XSSaggr, A, B, indexes, indexes_aggr, m_par, e_set, Copula, distrSS, compressionIndexes, priors,meas_error, meas_error_std)

Sample the posterior of the parameter vector using the Random-Walk Metropolis Hastings algorithm.

# Returns
- `draws::Array{Float64,2}`: `e_set.ndraws + e_set.burnin` sampled parameter vectors (row vectors)
- `posterior`: vector of posteriors for the respective draws
- `accept_rate`: acceptance rate
"""
function rwmh(xhat::Vector{Float64}, Σ::Symmetric{Float64,Array{Float64,2}}, n_par::NumericalParameters,
    Data::AbstractArray, D_miss::BitArray{2}, H_sel, XSSaggr, A, B, indexes, indexes_aggr, m_par, e_set, distrSS,
    compressionIndexes, priors, meas_error, meas_error_std)

    NormDist = MvNormal(zeros(length(xhat)), Σ)
    accept = 0
	accept_rate = 0.0
    draws = Matrix{Float64}(undef, e_set.ndraws + e_set.burnin, length(xhat))
    posterior = Vector{Float64}(undef, e_set.ndraws + e_set.burnin)
    draws[1, :] = xhat
    old_posterior, alarm = likeli(xhat, Data, D_miss, H_sel, XSSaggr, A, B, indexes, indexes_aggr,
                        m_par, n_par, e_set, distrSS, compressionIndexes, priors,
                        meas_error, meas_error_std)[3:4]
    posterior[1] = copy(old_posterior)
    proposal_draws = e_set.mhscale .* rand(NormDist, e_set.ndraws + e_set.burnin)
    for i = 2:e_set.ndraws + e_set.burnin
        xhatstar = draws[i-1, :] .+ proposal_draws[:, i]
        new_posterior, alarm, State2Control = likeli(xhatstar, Data, D_miss, H_sel, XSSaggr, A, B,
        indexes,indexes_aggr, m_par, n_par, e_set, distrSS,
        compressionIndexes, priors, meas_error, meas_error_std)[3:5]

        accprob = min(exp(new_posterior - old_posterior), 1.0)
        if alarm == false && rand() .<= accprob
            draws[i, :] = xhatstar
            posterior[i] = copy(old_posterior)
            old_posterior = new_posterior
            # @set! n_par.State2Control_save = State2Control
            accept += 1
        else
            draws[i, :] = draws[i-1, :]
            posterior[i] = posterior[i-1]
        end
        if mod(i, 200) == 0 || i == e_set.ndraws+ e_set.burnin
            print("-----------------------\n")
            print("Acceptance Rate: ", accept / i, "\n")
            print("Number of draws:", i,"\n")
            print("Parameters\n")
            print(draws[i, :],"\n")
            print("Posterior Likelihood:", old_posterior,"\n")
            print("-----------------------\n")
			accept_rate::Float64 = accept / i

        end

    end

	return draws, posterior, accept_rate
end

# draw overdispersed initial values for multi-chain RWMH
function multi_chain_init(xhat::Vector{Float64}, Σ::Symmetric{Float64,Array{Float64,2}}, n_par::NumericalParameters,
    Data::AbstractArray, D_miss::BitArray{2}, H_sel, XSSaggr,A,B, indexes, indexes_aggr,
	m_par, e_set, distrSS, compressionIndexes, priors, meas_error, meas_error_std)

	init_scale = 2 * e_set.mhscale # overdispersed initial values
	NormDist = MvNormal(zeros(length(xhat)), Σ)
	init_draw = Vector{Float64}(undef, length(xhat))
	init_success = false
	init_iter = 1
	while init_success == false && init_iter <= 100
		init_draw .= init_scale^2.0 .* rand(NormDist) .+ xhat

		alarm = likeli(init_draw, Data, D_miss, H_sel, XSSaggr, A, B, indexes, indexes_aggr,
                             m_par, n_par, e_set, distrSS, compressionIndexes, priors,
                             meas_error, meas_error_std)[4]
		if alarm == false
			init_success = true
		else
			init_iter += 1
		end
	end

	return init_draw, init_success, init_iter
end

function marginal_density(draws, posterior)
    ndraws, npars = size(draws)
	posterior_mode = maximum(posterior)
	d = Chisq(npars)
    θ_hat = mean(draws, dims = 1)[:]
    V_hat = cov(draws)
	inv_V_hat = inv(V_hat)

	marg_likeli_save = zeros(9)
	τ_iter = 1
	for τ = 0.1:0.1:0.9
		thresh = quantile(d, τ)
		const_terms = - 0.5 * npars * log(2 * pi) - 0.5 * logdet(V_hat) - log(τ)

		tmp = 0.0
	    for i = 1:ndraws
	        θ_dist = (draws[i, :] .- θ_hat)' * inv_V_hat * (draws[i, :] .- θ_hat)
	        if θ_dist <= thresh
	            log_f_θ = const_terms - 0.5 * θ_dist
	            tmp += exp(log_f_θ - posterior[i] + posterior_mode)
	        end
	    end
		marg_likeli_save[τ_iter] = posterior_mode - log(tmp / ndraws )
		τ_iter += 1
	end
    marg_likeli = mean(marg_likeli_save)

	return marg_likeli
end