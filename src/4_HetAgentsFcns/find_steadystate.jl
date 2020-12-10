@doc raw"""
    find_steadystate(m_par)

Find the stationary equilibrium capital stock.

# Returns
- `KSS`: steady-state capital stock
- `VmSS`, `VkSS`: marginal value functions
- `distrSS::Array{Float64,3}`: steady-state distribution of idiosyncratic states, computed by [`Ksupply()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`
"""
function find_steadystate(m_par)

BLAS.set_num_threads(Threads.nthreads())
# -------------------------------------------------------------------------------
## STEP 1: Find the stationary equilibrium for coarse grid
# -------------------------------------------------------------------------------
#-------------------------------------------------------
# Income Process and Income Grids
#-------------------------------------------------------
# Read out numerical parameters for starting guess solution with reduced income grid.
ny_temp                 = 5                             # Coarse initial income grid to obtain a first estimate of steady state capital stock 
grid_y, Π, bounds       = Tauchen(m_par.ρ_h,ny_temp)    # Income grid and transitions

# Include entrepreneurs into the income transitions
Π                       = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny_temp);
                            m_par.ι ./ ny_temp * ones(1, ny_temp) 1.0 .- m_par.ι]
grid_y                  = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2));
                            (m_par.ζ .+ m_par.ι)/m_par.ζ]

# Calculate expected level of human capital
Paux                    = Π^1000
H                       = (Paux[1,1:end-1]' * grid_y[1:end-1])

# Numerical parameters
n_par                   = NumericalParameters(ny = ny_temp+1, nm = 40, nk = 41, bounds_y = bounds, 
                                        grid_y = grid_y, Π = Π, H = H)
if n_par.verbose
    println("Finding equilibrium capital stock for coarse income grid")
end
n_par.grid_m[
  sum(n_par.grid_m.<0)] = 0.0 # Make sure zero is on the m grid (liquid asset)
n_par.mesh_m           .= repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)), outer=[1, n_par.nk, n_par.ny])

# Capital stock guesses
Kmax                    = 1.75 * ((m_par.δ_0 - 0.0025 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(1.0 / (m_par.α - 1.0))
Kmin                    = 1.0 * ((m_par.δ_0 - 0.0005 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(0.5 / (m_par.α - 1.0))

# a.) Define excess demand function
d(  K, 
    initial::Bool=true, 
    Vm_guess = zeros(1,1,1), 
    Vk_guess = zeros(1,1,1), 
    distr_guess = n_par.dist_guess
    )                   = Kdiff(K, n_par, m_par, initial, Vm_guess, Vk_guess, distr_guess)

# b.) Find equilibrium capital stock (multigrid on y,m,k)
KSS                     = CustomBrent(d, Kmin, Kmax)[1]
if n_par.verbose
    println("Capital stock is")
    println(KSS)
end
# -------------------------------------------------------------------------------
## STEP 2: Find the stationary equilibrium for final grid
# -------------------------------------------------------------------------------
if n_par.verbose
    println("Finding equilibrium capital stock for final income grid")
end
ny                      = n_par.ny_refined                      # use the final grid size for income from here on  
grid_y, Π, bounds       = Tauchen(m_par.ρ_h,ny)                 # Transition matrix, grid and integration bounds for income
grid_y                  = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2)); 
                            (m_par.ζ .+m_par.ι)/m_par.ζ]        # Add entrepreneur state

Π                       = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny);
                           m_par.ι ./ ny * ones(1,ny) 1.0 .- m_par.ι] # Add entrepreneur state

Paux                    = Π^1000
distr_y                 = Paux[1, :]                            # stationary income distribution
H                       = (distr_y[1:end-1]'*grid_y[1:end-1])   # average human capital

naggrstates             = length(state_names)                   # number of aggregate states
naggrcontrols           = length(control_names)                 # number of aggregate controls
naggr                   = length(aggr_names)                    # number of all aggregates

# Write changed parameter values to n_par
n_par                   = NumericalParameters(ny = ny+1, bounds_y = bounds, grid_y = grid_y, 
                        Π = Π, H = H, naggrstates = naggrstates, naggrcontrols = naggrcontrols,
                        aggr_names  = aggr_names, naggr = naggr,ϵ = 1e-10, 
                        )
n_par.grid_m[
        sum(n_par.grid_m.<0)
            ]           = 0.0
n_par.mesh_m           .= repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)),
                                    outer=[1, n_par.nk, n_par.ny])

# Find stationary equilibrium for refined economy
BrentOut                = CustomBrent(d, KSS*.95, KSS*1.05;tol = n_par.ϵ)
KSS                     = BrentOut[1]
VmSS                    = BrentOut[3][2]
VkSS                    = BrentOut[3][3]
distrSS                 = BrentOut[3][4]
if n_par.verbose
    println("Capital stock is")
    println(KSS)
end
return KSS, VmSS, VkSS, distrSS, n_par, m_par

end

