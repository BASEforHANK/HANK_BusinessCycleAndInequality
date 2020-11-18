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

# Define basic price functions
int_com(x)   = interest(x, 1 / m_par.μ, employment(x, 1 /
                        (m_par.μ*m_par.μw), m_par), m_par)
wage_com(x)  = wage(x, 1 / m_par.μ, employment(x, 1 /
               (m_par.μ*m_par.μw), m_par), m_par) * employment(x, 1 / m_par.μ, m_par)

profits(x)   = (1.0 - 1.0 ./ m_par.μ) * output(x, 1.0,
                employment(x, 1.0 ./ (m_par.μ*m_par.μw), m_par), m_par)

# -------------------------------------------------------------------------------
## STEP 1: Find the stationary equilibrium for coarse grid
# -------------------------------------------------------------------------------
#-------------------------------------------------------
# Income Process and Income Grids
#-------------------------------------------------------
# Read out numerical parameters for starting guess solution with reduced income grid.
ny_temp           = 7; # eigs in Ksupply quickly increases in runtime in ny (more than ny^2).
grid_y, Π, bounds = Tauchen(m_par.ρ_h,ny_temp) # Income grid and transitions
# Include entrepreneurs into the income transitions
Π               = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny_temp);
                    m_par.ι ./ ny_temp * ones(1, ny_temp) 1.0 .- m_par.ι]
grid_y          = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2));
                    (m_par.ζ .+ m_par.ι)/m_par.ζ]
# Calculate expected level of human capital
Paux            = Π^1000
H               = (Paux[1,1:end-1]' * grid_y[1:end-1])
# Numerical parameters
n_par           = NumericalParameters(ny = ny_temp+1, nm = 40, nk = 41, bounds_y = bounds, grid_y = grid_y, Π = Π, H = H)

n_par.grid_m[sum(n_par.grid_m.<0)] = 0.0 # Make sure zero is on the m grid (liquid asset)
n_par.mesh_m   .= repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)), outer=[1, n_par.nk, n_par.ny])

# Capital stock guesses
Kmax      = 1.75 * ((m_par.δ_0 - 0.0025 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(1.0 / (m_par.α - 1.0))
Kmin      = 1.0 * ((m_par.δ_0 - 0.0005 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(0.5 / (m_par.α - 1.0))
K         = range(Kmin, stop = Kmax, length = 8)
# a.) Define excess demand function
d(K, initial::Bool=true, Vm_guess = zeros(1,1,1), Vk_guess = zeros(1,1,1), distr_guess = n_par.dist_guess) =
        Kdiff(K, n_par, m_par, initial, Vm_guess, Vk_guess, distr_guess)

# b.) Find equilibrium capital stock (multigrid on y,m,k)
BrentOut  = CustomBrent(d, Kmin, Kmax)
KSS       = BrentOut[1]
distrSS   = BrentOut[3][4]
# -------------------------------------------------------------------------------
## STEP 2: Find the stationary equilibrium for final grid
# -------------------------------------------------------------------------------
# a) grid refinement
grid_y_old      = n_par.grid_y
grid_m_old      = n_par.grid_m
grid_k_old      = n_par.grid_k

ny              = n_par.ny_refined; # eigs in Ksupply quickly increases in runtime in ny
grid_y, Π, bounds= Tauchen(m_par.ρ_h,ny)
grid_y          = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2)); 
                        (m_par.ζ .+m_par.ι)/m_par.ζ]
Π               = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny);
                        m_par.ι ./ ny * ones(1,ny) 1.0 .- m_par.ι]
Paux            = Π^1000
distr_y         = Paux[1, :]
H               = (distr_y[1:end-1]'*grid_y[1:end-1]) 

naggrstates     = length(state_names)
naggr           = length(aggr_names)
naggrcontrols   = length(control_names)
# Write changed parameter values to n_par
n_par           = NumericalParameters(ny = ny+1, bounds_y = bounds, grid_y = grid_y, 
                        Π = Π, H = H, naggrstates = naggrstates, naggrcontrols = naggrcontrols,
                        aggr_names  = aggr_names, naggr = naggr,ϵ = 1e-10, 
                        )
n_par.grid_m[sum(n_par.grid_m.<0)] = 0.0
# Interpolate distribution on refined wealth-income grid
refined_dist    = mylinearinterpolate3(grid_m_old,grid_k_old, grid_y_old,
                                        distrSS, n_par.grid_m, n_par.grid_k, n_par.grid_y)
refined_dist    = refined_dist ./ sum(refined_dist, dims=(1, 2, 3))

# Add distribution guess
n_par           = NumericalParameters(ny = ny+1, bounds_y = bounds, grid_y = grid_y, 
                        Π = Π, H = H, naggrstates = naggrstates, naggrcontrols = naggrcontrols,
                        aggr_names  = aggr_names, naggr = naggr,ϵ = 1e-10, 
                        dist_guess  = refined_dist)
n_par.grid_m[sum(n_par.grid_m.<0)] = 0.0 
n_par.mesh_m   .= repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)),
                                    outer=[1, n_par.nk, n_par.ny])             
BrentOut                = CustomBrent(d, KSS*.9, KSS*1.2)
KSS                     = BrentOut[1]
VmSS                    = BrentOut[3][2]
VkSS                    = BrentOut[3][3]
distrSS                 = BrentOut[3][4]
return KSS, VmSS, VkSS, distrSS, n_par, m_par
end

