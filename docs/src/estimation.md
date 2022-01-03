# Estimation of parameters
!!! note
    Most of the code of this section is in the folder `6_Estimation`

## Settings
```@docs
EstimationSettings
```
The estimation settings (globally instantiated as `e_set` in `HANKEstim.jl`) manage
the following areas of the estimation:
- **the match between data and model**: `data_file` is a path to a .csv-file
    that contains quarterly observations of several variables (in columns), named in
    `observed_vars_input`. Missing data is denoted by `NaN`. If some column-names do not align with the model
    variables, `data_rename` is used. Level variables that should correspond to growth
    rates in the data can be selected in `growth_rate_select`. Measurement errors
    will be added for observables in `meas_error_input`
- **estimation of variances**: `shock_names` contain the aggregate shocks in the model,
    whose variances are estimated. `me_treatment` defines how measurement errors
    are treated: for `:fixed`, their variances are fixed by the data-variance, otherwise
    they are estimated either with `:bounded` uniform priors, or `:unbounded` priors 
    (see [`measurement_error()`](@ref)). For the latter case, the priors are set in
    `meas_error_distr`
- **numerical parameters**: the maximum number of iterations to find the mode of the
    likelihood (see [`mode_finding()`](@ref)) is set in `max_iter_mode`. `ndraws`, `burnin`
    and `mhscale` are parameters for the Random-Walk Metropolis Hastings algorithm (see [`rwmh()`](@ref))
- **estimation flags**: whether to estimate the model is set in `estimate_model`. `compute_hessian` determines whether the Hessian is computed after mode finding or set to an identity matrix (see [`mode_finding()`](@ref)). `multi_chain_init` sets whether multiple chains in the RWMH (see [`rwmh()`](@ref)) are started from an overdispersed posterior mode. All flags are set to `false` by default.

## Mode finding
```@docs
mode_finding
```
The main computations are the construction of the likelihood of the model parameters
and its maximization. We get the model parameters that are to be estimated,
together with their priors, from `m_par` (in addition to measurement error variances,
see [Settings](@ref)).

### The likelihood function
The function [`likeli()`](@ref) computes the log-likelihood of the model parameters `par`
in the following steps:

1. call [`SGU_estim()`](@ref) to derive the linear state-space representation of the model given `par`.
    Differently from [`SGU()`](@ref), differentiate only the system of *aggregate* equilibrium
    conditions with respect to *aggregate* variables, i.e. [`Fsys_agg()`](@ref). This is sufficient,
    as the estimated parameters do not enter in the heterogeneous agent part of the equilibrium system [^BBL].
    Then, update the derivatives `A` and `B` of the full model for aggregate variables and conditions,
    and compute the observation and state transition equations as in [`SGU()`](@ref)

2. delete rows of the observation equation that correspond to unobserved controls
    (the selector matrix `H_sel` is constructed in [`mode_finding()`](@ref)). Then, feed
    the linear state-space system, the data, and the variances of the structural and
    measurement shocks, into the Kalman filter (see [`kalman_filter()`](@ref)), which computes
    the log-likelihood

We find the maximizer of the likelihood function, as well as its Hessian at the maximum,
with the package `Optim`. Note that in order to obtain the Hessian, you need to set `e_set.compute_hessian = true`.

### Called functions
```@docs
likeli
SGU_estim
kalman_filter
measurement_error
```
## Bayesian estimation
```@docs
montecarlo
```
We use a Monte Carlo Markov Chain method, specifically the Random-Walk Metropolis Hastings ([`rwmh()`](@ref)) algorithm, to sample from the posterior probability distribution of the parameter vector. To obtain the
posterior likelihood of each draw, we call [`likeli()`](@ref), which evaluates the priors at `par` ([`prioreval()`](@ref)) and returns the log-posterior as a sum of the log-prior and the log-likelihood.

Given the draws from the posterior, we can analyze the probabilities of the parameters using the package `MCMCChains`. We take the average over the draws as our Bayesian estimate of the parameter vector, `par_final`. To obtain an estimate of the underlying state over the data sample period, we call [`likeli()`](@ref) with `par_final` and keyword `smoother=true` (this calls the Kalman smoother [`kalman_filter_smoother()`](@ref)). The result is stored in `smoother_output`, and saved with the other results in `e_set.save_posterior_file`.

### Called functions
```@docs
rwmh
prioreval
kalman_filter_smoother
```

[^BBL]:
    For details, see the paper [Shocks, Frictions, and Inequality in US Business Cycles](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=14364)