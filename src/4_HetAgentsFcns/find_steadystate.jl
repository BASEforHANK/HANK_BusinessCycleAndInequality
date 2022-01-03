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

# -------------------------------------------------------------------------------
## STEP 1: Find the stationary equilibrium for coarse grid
# -------------------------------------------------------------------------------
#-------------------------------------------------------
# Income Process and Income Grids
#-------------------------------------------------------
# Read out numerical parameters for starting guess solution with reduced income grid.

# Numerical parameters 
n_par   = NumericalParameters(m_par = m_par, ny = 4, nm = 10, nk = 10,  ϵ = 1e-6)
if n_par.verbose
    println("Finding equilibrium capital stock for coarse income grid")
end
# Capital stock guesses
Kmax   = 1.75 * ((m_par.δ_0 - 0.0025 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(1.0 / (m_par.α - 1.0))
Kmin   = 0.5 * ((m_par.δ_0 - 0.0025 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(1.0 / (m_par.α - 1.0))

# a.) Define excess demand function
d(  K, 
    initial::Bool=true, 
    Vm_guess = zeros(1,1,1), 
    Vk_guess = zeros(1,1,1), 
    distr_guess = n_par.dist_guess) = Kdiff(K, n_par, m_par, initial, Vm_guess, Vk_guess, distr_guess)

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
# Write changed parameter values to n_par
n_par                   = NumericalParameters(m_par = m_par, naggrstates = length(state_names), naggrcontrols = length(control_names),
                                              aggr_names  = aggr_names, distr_names = distr_names)

# Find stationary equilibrium for refined economy
BrentOut                = CustomBrent(d, KSS*.8, KSS*1.2;tol = n_par.ϵ)
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

