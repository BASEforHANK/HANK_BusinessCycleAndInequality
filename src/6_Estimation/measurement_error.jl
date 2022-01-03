@doc raw"""
    measurement_error(Data,observed_vars,e_set)

Build measurement error.

# Arguments
- `Data`: matrix of observables [nobs * nvar]
- `observed_vars`: vector of observed variable names [nvar * 1]
- `e_set::EstimationSettings`

# Returns
- `meas_error`: ordered dictionary of measurement errors linked to observables
- `meas_error_prior`: corresponding priors for measurement errors
- `meas_error_std`: standard deviations of observables with measurement error
"""
function measurement_error(Data, observed_vars, e_set)

    if !isempty(e_set.meas_error_input)

        # find correct positions for measurement error
        meas_error_index = Vector{Int}(undef,length(e_set.meas_error_input))
        for i in eachindex(e_set.meas_error_input)
            meas_error_index[i] = findall(x->x==e_set.meas_error_input[i], observed_vars)[1]
        end
        meas_error = OrderedDict(zip(e_set.meas_error_input, meas_error_index))
        
        # create measurement error according to selected treatment
        if e_set.me_treatment == :unbounded
            # inverse gamma prior
            meas_error_prior = e_set.meas_error_distr
            meas_error_std = e_set.me_std_cutoff * ones(length(e_set.meas_error_input))
        elseif e_set.me_treatment == :bounded || e_set.me_treatment == :fixed
            # data dependent hard upper bound on measurement error standard deviation
            meas_error_prior = Array{Uniform{Float64}}(undef,length(e_set.meas_error_input))
            meas_error_std = Vector{Float64}(undef, length(e_set.meas_error_input))
            m_iter = 1
            for (k, v) in meas_error # read out position of measurement errors
                meas_error_std[m_iter] = e_set.me_std_cutoff * std(skipmissing(Data[:, v]))
                meas_error_prior[m_iter] = Uniform(0.0, meas_error_std[m_iter])
                m_iter += 1
            end
        else
            error("ME treatment not implemented")
        end
    else
        # in case of no measurement error
        meas_error = OrderedDict{Symbol, Int}()
        meas_error_prior = repeat([InverseGamma(ig_pars(0.0005, 0.001^2)...)], 0)
        meas_error_std = Vector{Float64}(undef, 0)
    end

    return meas_error, meas_error_prior, meas_error_std
end
