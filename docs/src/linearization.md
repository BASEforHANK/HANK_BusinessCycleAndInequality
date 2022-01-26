# Linear perturbation around steady state

!!! note
    The main functions of this section are in the folder `5_LinearizationFunctions`.

The model is linearized with respect to aggregate variables. For this,
we write the equilibrium conditions in the form of
``F(X,X')=0``, where ``X`` and ``X'`` are (expected) deviations from steady state
in two successive periods. Applying the total differential yields
``A*X' = - B*X``, where ``A``,``B`` are the first derivatives of ``F`` with respect
to ``X'``,``X``. In the standard setting, we use the generalized Schur decomposition [^Klein]
to transform this equation into a linearized observation equation ``d = gx*k`` and
a linearized state transition equation ``k' = hx*k``, where ``k`` is a vector of the
*state* variables and ``d`` is a vector of the *control* variables (``X = \begin{bmatrix} k \\ d \end{bmatrix}``).

In our code, ``F`` is implemented as [`HANKEstim.Fsys()`](@ref), while differentiating and
solving for ``gx`` and ``hx`` is done in [`HANKEstim.SGU()`](@ref), and [`linearize_full_model()`](@ref)
returns the results as a `struct` `LinearResults`:
```@docs
linearize_full_model
```
## Overview of `SGU()`
```@docs
HANKEstim.SGU
```
The function executes the following steps:

- generate devices to retrieve distribution and marginal value functions from
    compressed states/controls (`Γ` and `DC`,`IDC`)
- calculate the first derivative of [`HANKEstim.Fsys()`](@ref) with respect to `X` and `XPrime`.
    We use automatic differentiation (implemented in Julia by the package `ForwardDiff`).
    Partial derivatives are calculated using the `ForwardDiff.jacobian()` function.
    We exploit that some partial derivatives have known values (contemporaneous marginal value
    functions and the future marginal distributions) and set them directly instead of calculating them [^BL].

- compute linear observation and state transition equations using the [`HANKEstim.SolveDiffEq()`](@ref) function

## Overview of `SolveDiffEq()'
```@docs
HANKEstim.SolveDiffEq
```
- compute linear observation and state transition equations. The solution algorithm is set
    in `n_par.sol_algo`, with the options `:schur` (mentioned above) and `:litx` [^lit]. The results are matrices that map contemporaneous states to controls [`gx`],
    or contemporaneous states to future states [`hx`]


## Overview of `Fsys()`
```@docs
HANKEstim.Fsys
```
The function [`HANKEstim.Fsys()`](@ref) proceeds in the following way:
1. set up vector `F`, that contains the errors to all equilibrium conditions. There are as many conditions
    as deviations from steady state (length of `X`,`XPrime`), and conditions are indexed with
    respective model variable in `IndexStruct` `indexes`
2. generate locally all aggregate variables (for both periods) using [`@generate_equations`](@ref)
3. construct the full-grid marginal distributions, marginal value functions, and the copula
    from the steady-state values and the (compressed) deviations (for the copula, the selection of DCT
    coefficients that can be perturbed ensures that also the perturbed function is a copula)
4. write all equilibrium condition-errors with respect to *aggregate* variables to `F`, using
    [`HANKEstim.Fsys_agg()`](@ref)
5. compute optimal policies with [`HANKEstim.EGM_policyupdate()`](@ref), given
    future marginal value functions, prices, and individual incomes. Infer present marginal
    value functions from them (envelope theorem) and set the difference to assumed present
    marginal value functions (in terms of their compressed deviation from steady state)
    as equilibrium condition-errors (*backward iteration of the value function*)
6. compute future marginal distributions and the copula (on the copula grid) from previous distribution and optimal asset policies. Interpolate when necessary. Set difference to assumed future marginal distributions and copula values on the copula nodes as equilibrium condition-errors (*forward iteration of the distribution*)
7. compute distribution summary statistics with [`HANKEstim.distrSummaries()`](@ref) and write
    equilibrium conditions with their respective (control) variables
8. return `F`

Note that the copula is treated as the sum of two interpolants. An interpolant based on the steady-state distribution using the full steady-state marginals as a grid and a "deviations"-function that is defined on the copula grid generated in `prepare_linearization()`. The actual interpolation is carried out with [`HANKEstim.myinterpolate3()`](@ref). Default setting is trilinear interpolation, the code also allows for 3d-Akima interpolation.

### Called functions / macros
```@docs
@generate_equations
HANKEstim.Fsys_agg
HANKEstim.myinterpolate3
```

[^Klein]:
    See the paper [Using the generalized Schur form to solve a multivariate linear rational expectations model](https://www.sciencedirect.com/science/article/pii/S0165188999000457) by Paul Klein (JEDC 2000)

[^BL]:
    Contemporaneous marginal value functions are irrelevant for optimal decisions, so
    its effect on other model variables is 0. Due to a rich enough set of prices, the future distribution
    directly only affects the Fokker-Planck equation. For details, see the paper
    [Solving heterogeneous agent models in discrete time with many idiosyncratic states by perturbation methods](https://cepr.org/active/publications/discussion_papers/dp.php?dpno=13071#)

[^lit]:
    Invoking the Implicit Function Theorem, there exist functions ``g`` and ``h`` such that
    ``F\left(\begin{pmatrix} k \\ g(k) \end{pmatrix},\begin{pmatrix} h(k) \\ g(h(k)) \end{pmatrix}\right)=0``.
    Totally differentiating by ``k`` yields ``B \begin{pmatrix}\mathbb{I}\\ Dg \end{pmatrix}+A \begin{pmatrix}\mathbb{I}\\ Dg \end{pmatrix} Dh = 0``. The `:lit`-algorithm solves this equation for ``Dg`` and ``Dh`` iteratively.
