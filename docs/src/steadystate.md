# Computation of the steady state and dimensionality reduction
!!! note
    Most of the code of this section is in the folder `4_HetAgentsFcns`, except for  `prepare_linearization()` which is in `5_LinearizationFunctions`.

The model features uninsured income shocks ``y`` (by assumption, all workers supply the same
efficiency units of labor [^BBL], so idiosyncratic productivity shocks translate
to income shocks) and two assets, bonds ``m`` and illiquid capital ``k``. Entrepreneurs
(last income-state) receive no labor income, but firm profits, while workers additionally
receive labor union profits.                                                                                                

The steady state equilibrium contains marginal value functions ``V_m`` and ``V_k``
on a three-dimensional grid ``(m \times k \times y)`` and the ergodic joint distribution
over these idiosyncratic states. We do dimensionality reduction [^BL] by applying
the Discrete Cosine Transformation to the marginal value functions and approximating
the joint distribution with a copula and state-dependent marginals.

The main functions are [`HANKEstim.find_steadystate()`](@ref) and [`HANKEstim.prepare_linearization()`](@ref):

## Overview of `find_steadystate`
```@docs
HANKEstim.find_steadystate
```
The function takes the parameter `struct` `ModelParameters` as input `m_par` (see [Parameters](@ref)).

To find the stationary equilibrium, we proceed in roughly the following steps:

1. instantiate the parameter `struct` `NumericalParameters` as `n_par` (see [Parameters](@ref)).
   Within the struct, we set the number of income states [`ny`] and use the [`HANKEstim.Tauchen()`](@ref) method to obtain a grid and a transition matrix of income, given the autocorrelation of the income process [`m_par.ρ_h`]. Then, include entrepreneurial state.
2. find equilibrium capital stock (by finding a root of [`HANKEstim.Kdiff()`](@ref)), where
    the supply of capital by households is calculated in [`HANKEstim.Ksupply()`](@ref),
    which uses the Endogenous Grid Method (see [`HANKEstim.EGM_policyupdate`](@ref))
    to iteratively obtain optimal policies and marginal value functions

## Overview of `prepare_linearization`
```@docs
HANKEstim.prepare_linearization
```
We first calculate other equilibrium quantities and produce distributional summary statistics ([`HANKEstim.distrSummaries()`](@ref)). Next, we reduce the dimensionality:

1. compute coefficients of the Chebyshev polynomials that serve as basis functions
    for ``V_m`` and ``V_k``, using the Discrete Cosine Transformation (Julia-package
    `FFTW`), and retain those that explain the most of its variance, up to 
    `100*(1-n_par.reduc)` percent. Save their indices in `compressionIndexes`
2. prepare a node mesh on which the time-varying linear interpolant of the copula
     is defined. The grid in each ``m``, ``k``, and ``y`` dimension is selected 
     such that each resulting bin holds approximately the same share of the
     respective aggregate variable. 

Lastly, we collect the steady state values of all model variables in the 
vector `XSS` (see [`@writeXSS`](@ref)). The *state* variables consist of
the marginal distributions over ``m``, ``k`` and ``y`` and the aggregate state variables
(collected in `state_names`). The *control* variables consist of the steady state
marginal value functions (over the full grid) and the aggregate control variables
(collected in `control_names`; these vectors are defined in the main script `HANKEstim.jl`).

While the steady state marginal value functions have full dimensionality,
in the vectors that collect *deviations* from steady state (in [`HANKEstim.Fsys()`](@ref), those are `X` and `XPrime`)
only the coefficients of the most important Chebyshev polynomials are saved.
Additionally, the deviations of the marginal distributions are saved with one entry short of
the grid size, since the marginals are restricted to sum up to 1.
We manage this by creating the `struct` `indexes` (using [`@make_fn`](@ref)),
that has two fields for each variable: steady state value and deviation.

We also construct the vector `XSSaggr` and the `struct` `indexes_aggr`,
which are similar to the above but only store (and manage) aggregate variables.
This is useful for differentiating only with respect to aggregate variables
in the estimation part (see [`HANKEstim.SGU_estim()`](@ref)).

!!! warning
    Be sure that you edit `prepare_linearization()` and not `prepare_linearization_generated()` which will be overwritten by the model parser based on `prepare_linearization()`.

## Parameters
The model parameters for the steady state have to be calibrated. We set them
in the `struct` `ModelParameters`. It also contains all other parameters that
are estimated, including the stochastic process-parameters for the aggregate
shocks.
```@docs
ModelParameters
```
The numerical parameters contain the grid (and the meshes) on which the
stationary equilibrium is solved, discretization results of [`HANKEstim.find_steadystate()`](@ref) 
like the transition matrix of income and the joint distribution, and other
parameters that determine the numerical approximation or solution technique,
like `reduc` or `sol_algo`.
```@docs
NumericalParameters
```
## Find stationary equilibrium: functions
```@docs
HANKEstim.Tauchen
HANKEstim.Kdiff
HANKEstim.Ksupply
HANKEstim.EGM_policyupdate
HANKEstim.distrSummaries
```

## Collect variables: macros
```@docs
@writeXSS
@make_fn
@make_fnaggr
@make_struct
@make_struct_aggr
```

[^BBL]:
    For details, see the paper [Shocks, Frictions, and Inequality in US Business Cycles](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=14364)
[^BL]:
    For details, see the paper
    [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=13071#)