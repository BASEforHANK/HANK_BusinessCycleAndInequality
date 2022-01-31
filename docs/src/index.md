# HANKEstim.jl Documentation
## Introduction
This manual documents the Julia module **HANKEstim**, that implements the solution and Bayesian likelihood
estimation of a heterogeneous-agent New-Keynesian (HANK) model. It accompanies the paper
[Shocks, Frictions, and Inequality in US Business Cycles](https://www.benjaminborn.de/publication/bbl_inequality_2020/).

!!! note
    The paper is currently under revision and the updated toolbox does not replicate the results in the linked version of the paper.

## First steps
The module runs with Julia 1.7.1. We recommend to use [Julia for VSCode IDE](https://www.julia-vscode.org) as a front-end to Julia. To get started with the toolbox, simply download or clone the folder, e.g. via `git clone`, `cd` to the project directory and call
```julia-repl
(v1.7) pkg> activate .

(HANK_BusinessCycleAndInequality) pkg> instantiate
```
This will install all needed packages. For more on Julia environments, see [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project).

!!! warning
    Before you activate the environment, make sure that you are in the main directory, in which the `Manifest.toml` and `Project.toml` files are located. In case you accidentally activated the environment in a subfolder, empty `.toml` files will be created that you need to delete before proceeding in the correct folder.


For an introduction, it is easiest to use the Julia script `script.jl` in the `src` folder. Make sure that the folder is the present working directory and that the bottom bar in VSCode shows `Julia env: HANK_BusinessCycleAndInequality`.[^1] At the top of the script file, we pre-process some user input regarding the aggregate model and the steady state (see below) and write them them into the respective functions in the folder `5_LinearizationFunctions\generated_fcns`. This has to be done before the HANKEstim module, defined in `HANKEstim.jl`, is loaded via
```
push!(LOAD_PATH, pwd())
using HANKEstim
```
```@meta
# !!! note
#
#    Instead of pushing the current directory to `LOAD_PATH` at runtime, one can also move the folder `HANKEstim` to
#    the place where packages are stored in the local Julia environment.
```
`HANKEstim.jl` is the key module file as it loads in the code base, sets up structures, and exports a number of functions and macros.

The provided `script.jl` then shows how a typical estimation proceeds in three main steps. First, we solve the steady state of the model, and reduce the dimensionality of the state space [^BL]. Secondly, we compute the linearized dynamics of the model around the steady state. Thirdly, we construct the likelihood of the model parameters given data, and use Bayesian methods to estimate them. More details on the three steps are provided in the menu on the left. `script.jl` also provides an example on how to plot some impulse response functions from the model.

### Setting up your model

To define the aggregate part of the model, include the aggregate model block in `1_Model\input_aggregate_model.jl`. The model variables are divided into *states* (distribution, productivity, ...) and
*controls* (consumption policy or marginal utilities, prices, aggregate capital stock, ...). The aggregate variables (i.e. excluding the distribution and marginal utilities) are defined
 in `1_Model\include_aggregate_names` and their steady states in `1_Model\input_aggregate_steady_state`. Include model parameters in `struct ModelParameters` in `1_Model\Parameters.jl`.

The file `Parameters.jl` contains three structures to provide model parameters, numerical parameters, and estimation settings. In additon, it contains two macros that automatically create structures that contain the model variables.

The model parameters for the steady state have to be calibrated. We set them in the `struct` [`ModelParameters`](@ref). It also contains all other parameters that are estimated, including the stochastic process-parameters for the aggregate shocks. Each model parameter has a line of code. It starts with the parameter name as it is used in the code and a default value. The next two entries are its ascii name and its name for LaTeX output. The fourth entry is the prior if the parameter is to be estimated. Please see the [Distributions.jl](https://github.com/JuliaStats/Distributions.jl)-package for available options. The fifth entry is a Boolean whether the parameter should be estimated (`true`) or not (`false`)


### Steady state and first dimensionality reduction
The command
```
sr_full = compute_steadystate(m_par)
```
calls the functions [`HANKEstim.find_steadystate()`](@ref) and [`HANKEstim.prepare_linearization()`](@ref) and saves their returns in an instance `sr_full` of the `struct` `SteadyResults`.
`sr_full` contains vectors of the steady-state variables (together with index-vectors to reference them by name) and
the steady-state distribution of income and assets. It also contains the compressed marginal value functions(using DCTs on the steady-state value functions).

!!! tip
    `sr_full` may be saved to the local file system by calling
    ```
    @save "7_Saves/steadystate.jld2" sr_full
    ```
    and can be loaded for a future session with
    ```
    @load "7_Saves/steadystate.jld2" sr_full
    ```

### Linearize full model
After computing the steady state and saving it in the `SteadyResults`-struct named `sr_full`,
```
lr_full = linearize_full_model(sr_full, m_par)
```
computes the linear dynamics of the model around the steady state (in the background, this calls [`HANKEstim.SGU_estim()`](@ref)) and saves a state-space representation
in the instance `lr_full` of the `struct` `LinearResults` (see [`linearize_full_model()`](@ref)).

Linearization of the full model takes a few seconds. The resulting state space is, because the copula and the value functions are treated fully flexible in this first step, relatively large. As a result, also computing the first-order dynamics of this model takes a few seconds as well.

### Model reduction
This large state-space representation can however be reduced substantially. For this purpose, run  
```
sr_reduc    = model_reduction(sr_full, lr_full, m_par)
```
which calculates the unconditional covariance matrix of all state and control variables and rewrites the coefficients of the value functions and the copula as linear combinations of some underlying factors. Only those factors that have eigenvalues above the precision predefined in `sr_full.n_par.compress_critC` and `sr_full.n_par.compress_critS` are retained.
!!! warning
    After model reduction, `sr_reduc.indexes_r` contains the indexes that map correctly into the states/controls used in `LOMstate` and `State2Control`.


### Model solution after a parameter change / after reduction
This smaller model (or any model after a parameter change that doesn't affect the steady state) can be solved quickly using a factorization result from [^BBL] running
```
lr_reduc    = update_model(sr_reduc, lr_full, m_par)
```
In the background, this calls [`HANKEstim.SGU_estim()`](@ref), which only updates the Jacobian entries that regard the **aggregate** model. (Note that both [`HANKEstim.SGU()`](@ref) and [`HANKEstim.SGU_estim()`](@ref) call [`HANKEstim.SolveDiffEq()`](@ref) to obtain a solution to the linearized difference equation.)

This model update step takes about 200ms on a standard computer for a medium size resolution.

### Estimation of model parameters
Having obtained `SteadyResults` `sr_reduc` and `LinearResults` `lr_reduc`, the command
```
er_mode = find_mode(sr_reduc, lr_reduc, m_par)
```
computes the mode of the likelihood, i.e., the parameter vector that maximizes the probability of
observing the data given the model, and saves the results in `er_mode`, an instance of `struct` `EstimResults`
(see [`HANKEstim.mode_finding()`](@ref)). We use the Kalman filter to compute the likelihood, and the package
`Optim` for optimization. Settings for the estimation can adjusted in the `struct` [`EstimationSettings`](@ref).

!!! warning
    By default, the flag `estimate_model` in the `struct` [`EstimationSettings`](@ref) is set to `false`. Depending on the computing power available, finding the mode of the likelihood can take several hours to run through. The mode finder might also seem frozen after finishing the optimization but the computation of the Hessian for the large model is involved and can take a long time for the large model. For instructional purposes, we therefore set `e_set.compute_hessian = false` by default and load the Hessian from a save file. For a proper estimation, this has to be set to true. We also save an intermediate step before computing the Hessian in case you are only interested in the mode itself.

Lastly,
```
montecarlo(sr_reduc, lr_reduc, er_mode, m_par)
```
uses a Monte Carlo Markov Chain method to trace out the posterior probabilites of the estimated parameters.
The final estimates (and further results) are saved in a file with the name given by the field `save_posterior_file`
in the `struct` `EstimationSettings` (instantiated in `e_set`).

!!! note
    The module `HANKEstim` creates the estimation settings `e_set` in its main script (when it is initialized),
    so changes to the `struct` `EstimationSettings` are only effective *before* `using HANKEstim`. Make sure
    that all file paths specified in [`EstimationSettings`](@ref) are correct relative to your script's position.
[^1]:
    If you use a different editor, make sure that the environment is correctly set, as otherwise the instantiated packages might not be found.
[^BBL]:
    See the paper [Shocks, Frictions, and Inequality in US Business Cycles](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=14364)
[^BL]:
    For a description of the solution methods applied here, see the paper
    [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=13071#)
