# HANKEstim.jl Documentation
## Introduction
This manual documents the Julia module **HANKEstim**, that implements a Bayesian likelihood
estimation of a heterogeneous-agent New-Keynesian (HANK) model. It accompanies the paper
[Shocks, Frictions, and Inequality in US Business Cycles](https://www.benjaminborn.de/publication/bbl_inequality_2020/).
## Basic usage
The module runs with Julia 1.4. To use the toolbox in a Julia script, move the script to the directory that contains the folder `HANKEstim`, make sure that the folder is the present working directory and write
```
push!(LOAD_PATH, pwd())
using HANKEstim
```
!!! note

    Instead of pushing the current directory to `LOAD_PATH` at runtime, one can also move the folder `HANKEstim` to
    the place where packages are stored in the local Julia environment.

You might also have to install a number of packages that are listed at the top of `HANKEstim.jl`.

The estimation proceeds in three main steps. First, we solve the steady state of the model, and reduce the
dimensionality of the state space [^BL]. Secondly, we compute the linearized dynamics of the model
around the steady state. Thirdly, we construct the likelihood of the model parameters given data,
and use Bayesian methods to estimate them.

### Setting up your model

To define the aggregate part of the model, include the aggregate model block in `input_aggregate_model.jl`. States and controls have to be listed in `include_aggregate_names` and their steady states in `input_aggregate_steady_state`. Include model parameters in `struct ModelParameters` in `Struct.jl`.

!!! note

    For replication purposes, currently the estimated parameters in `Struct.jl` are overwritten in `find_SS()` by the estimated values from the paper.



### Steady state and dimensionality reduction
The command
```
sr = compute_steadystate()
```
calls the function [`find_SS()`](@ref) and saves its returns in an instance `sr` of the `struct` `SteadyResults`.
`sr` contains vectors of the steady state variables (together with index-vectors to reference them by name),
the steady-state distribution of income and assets, and devices to retrieve the full states from the
compressed state vectors.

!!! tip
    `sr` may be saved to the local file system by calling
    ```
    save_steadystate(sr)
    ```
    and can be loaded for a future session with
    ```
    sr = load_steadystate()
    ```
    These functions are necessary as the functions `save()` and `load()` from the `JLD2` package
    show incompatibilities with the more complex data structures in `SteadyResults`, especially
    the `struct` `NumericalParameters`.

The model variables are divided into *states* (distribution, productivity, ...) and
*controls* (consumption policy or marginal utilities, prices, aggregate capital stock, ...).
The aggregate variables (i.e. excluding the distribution and marginal utilities) are defined
in the vectors `state_names` and `control_names` in the main script `HANKEstim.jl`.

### Linearize full model
After obtaining `SteadyResults` `sr`,
```
lr = linearize_full_model(sr)
```
computes the linear dynamics of the model around the steady state and saves a state-space representation
in the instance `lr` of the `struct` `LinearResults` (see [`linearize_full_model()`](@ref)).

While solving for the first order dynamics of the full model takes only a few seconds (after having been compiled, i.e. after the first run),
a factorization result [^BBL] makes it possible to only solve the *aggregate* part of the model
when estimating the parameters (see [`SGU_estim()`](@ref)), which significantly reduces its computation time.

### Estimation of model parameters
Having obtained `SteadyResults` `sr` and `LinearResults` `lr`, the command
```
er = find_mode(sr,lr)
```
computes the mode of the likelihood, i.e. the parameter vector that maximizes the probability of
observing the data given the model, and saves the results in `er`, an instance of `struct` `EstimResults`
(see [`mode_finding()`](@ref)). We use the Kalman filter to compute the likelihood, and the package
`Optim` for optimization.

!!! warning
    Depending on the computing power available, finding the mode of the likelihood can take several hours to run through.

Lastly,
```
montecarlo(sr,lr,er)
```
uses a Monte Carlo Markov Chain method to trace out the posterior probabilites of the estimated parameters.
The final estimates (and further results) are saved in a file with the name given by the field `save_posterior_file`
in the `struct` `EstimationSettings` (instantiated in `e_set`).

!!! note
    The module `HANKEstim` creates the estimation settings `e_set` in its main script (when it is initialized),
    so changes to the `struct` `EstimationSettings` are only effective *before* `using HANKEstim`. Make sure
    that all file paths specified in [`EstimationSettings`](@ref) are correct relative to your script's position.

[^BBL]:
    See the paper [Shocks, Frictions, and Inequality in US Business Cycles](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=14364)
[^BL]:
    For a description of the solution methods applied here, see the paper
    [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=13071#)
