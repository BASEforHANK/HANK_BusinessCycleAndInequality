# Computation of the steady state and dimensionality reduction
!!! note
    Most of the code of this section is in the folder `4_HetAgentsFcns`.

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

The main functions are [`find_steadystate()`](@ref) and [`prepare_linearization()`](@ref):

## Overview of `find_steadystate`
```@docs
find_steadystate
```
The function takes the parameter `struct` `ModelParameters` as input `m_par` (see [Parameters](@ref)).

To find the stationary equilibrium, we proceed in roughly the following steps:

1. instantiate the parameter `struct` `NumericalParameters` as `n_par` (see [Parameters](@ref)).
   Within the struct, we set the number of income states [`ny`] and use the [`Tauchen()`](@ref) method to obtain a grid and a transition matrix of income, given the autocorrelation of the income process [`m_par.ρ_h`]. Then, include entrepreneurial state.
2. find equilibrium capital stock (by finding a root of [`Kdiff()`](@ref)), where
    the supply of capital by households is calculated in [`Ksupply()`](@ref),
    which uses the Endogenous Grid Method (see [`EGM_policyupdate`](@ref))
    to iteratively obtain optimal policies and marginal value functions

## Overview of `prepare_linearization`
```@docs
prepare_linearization
```
We first calculate other equilibrium quantities and produce distributional summary statistics ([`distrSummaries()`](@ref)). Next, we reduce the dimensionality:

1. compute coefficients of the Chebyshev polynomials that serve as basis functions
    for ``V_m`` and ``V_k``, using the Discrete Cosine Transformation (Julia-package
    `FFTW`), and retain those that explain the most of its variance, up to 
    `100*(1-n_par.reduc)` percent. Save their indices in `compressionIndexes`
2. compute the Copula as a function that maps three marginal
    distributions to a linear interpolation of the joint distribution on its
    marginals (see [`mylinearinterpolate3()`](@ref))

Lastly, we collect the steady state values of all model variables in the 
vector `XSS` (see [`@writeXSS`](@ref)). The *state* variables consist of
the marginal distributions over ``m``, ``k`` and ``y`` and the aggregate state variables
(collected in `state_names`). The *control* variables consist of the steady state
marginal value functions (over the full grid) and the aggregate control variables
(collected in `control_names`; these vectors are defined in the main script `HANKEstim.jl`).

While the steady state marginal value functions have full dimensionality,
in the vectors that collect *deviations* from steady state (in [`Fsys()`](@ref), those are `X` and `XPrime`)
only the coefficients of the most important Chebyshev polynomials are saved.
Additionally, the deviations of the marginal distributions are saved with one entry short of
the grid size, since the marginals are restricted to sum up to 1.
We manage this by creating the `struct` `indexes` (using [`@make_fn`](@ref)),
that has two fields for each variable: steady state value and deviation.

We also construct the vector `XSSaggr` and the `struct` `indexes_aggr`,
which are similar to the above but only store (and manage) aggregate variables.
This is useful for differentiating only with respect to aggregate variables
in the estimation part (see [`SGU_estim()`](@ref)).

## Parameters
The model parameters for the steady state have to be calibrated. We set them
in the `struct` `ModelParameters`. It also contains all other parameters that
are estimated, including the stochastic process-parameters for the aggregate
shocks.
```@docs
ModelParameters
```
The numerical parameters contain the grid (and the meshes) on which the
stationary equilibrium is solved, discretization results of [`find_steadystate()`](@ref) 
like the transition matrix of income and the joint distribution, and other
parameters that determine the numerical approximation or solution technique,
like `reduc` or `sol_algo`.
```@docs
NumericalParameters
```
## Find stationary equilibrium: functions
```@docs
Tauchen
Kdiff
Ksupply
EGM_policyupdate
distrSummaries
```
## Dimensionality reduction: functions
```@docs
mylinearinterpolate3
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