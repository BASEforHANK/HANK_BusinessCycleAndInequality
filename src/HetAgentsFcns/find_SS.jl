@doc raw"""
    find_SS(refined_ny)

Find the stationary equilibrium using a final resolution of `refined_ny` for income.

# Returns
- `XSS::Array{Float64,1}`, `XSSaggr::Array{Float64,1}`: steady state vectors produced by [`@writeXSS()`](@ref)
- `indexes`, `indexes_aggr`: `struct`s for accessing `XSS`,`XSSaggr` by variable names, produced by [`@make_fn()`](@ref),
        [`@make_fnaggr()`](@ref)
- `compressionIndexes::Array{Array{Int,1},1}`: indexes for compressed marginal value functions (``V_m`` and ``V_k``)
- `Copula(x,y,z)`: function that maps marginals `x`,`y`,`z` to approximated joint distribution, produced by
        [`mylinearinterpolate3()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`
- `distrSS::Array{Float64,3}`: steady state distribution of idiosyncratic states, computed by [`Ksupply()`](@ref)
- `CDF_SS`, `CDF_m`, `CDF_k`, `CDF_y`: cumulative distribution functions (joint and marginals)
"""
function find_SS()

BLAS.set_num_threads(Threads.nthreads())

# global m_par, n_par, CDF_m, CDF_k, CDF_y
# load estimated parameter set
m_par           = ModelParameters( )
@load e_set.mode_start_file par_final parnames
par = par_final[1:length(parnames)]
if e_set.me_treatment != :fixed
  m_par = Flatten.reconstruct(m_par, par[1:length(par) - length(e_set.meas_error_input)])
else
  m_par = Flatten.reconstruct(m_par, par)
end

# Define basic price functions
int_com(x)      = interest(x, 1 / m_par.μ, employment(x, 1 /
                        (m_par.μ*m_par.μw), m_par), m_par)
wage_com(x)     = wage(x, 1 / m_par.μ, employment(x, 1 /
                 (m_par.μ*m_par.μw), m_par), m_par) *
                  employment(x, 1 / m_par.μ, m_par)

profits(x)      = (1.0 - 1.0 ./ m_par.μ) * output(x, 1.0,
                    employment(x, 1.0 ./
                    (m_par.μ*m_par.μw), m_par), m_par)

# Read out numerical parameters for starting guess solution with reduced income grid.
ny              = 5; # eigs in Ksupply quickly increases in runtime in ny (more than ny^2).
grid_y, Π, bounds = Tauchen(m_par.ρ_h,ny) # Income grid and transitions
# Include entrepreneurs into the income transitions
Π               = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny);
                    m_par.ι ./ ny * ones(1,ny) 1.0 .- m_par.ι]
grid_y          = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2));
                    (m_par.ζ .+ m_par.ι)/m_par.ζ]
# Calculate expected level of human capital
Paux            = Π^1000
H               = (Paux[1,1:end-1]'*grid_y[1:end-1])
# Numerical parameters
n_par           = NumericalParameters(ny = ny+1, kmax = 1500, bounds_y = bounds,
                    mmin =-6.6, mmax = 1000, ϵ = 1e-5, grid_y = grid_y, Π=Π,H=H)

n_par.grid_m[sum(n_par.grid_m.<0)] = 0.0 # Make sure zero is on the m grid (liquid asset)
n_par.mesh_m   .= repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)),
                                    outer=[1, n_par.nk, n_par.ny])

# -------------------------------------------------------------------------------
## STEP 1: Find the stationary equilibrium
# -------------------------------------------------------------------------------
# Capital stock guesses
Kmax      = 2.0 * ((m_par.δ_0 - 0.0025 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(1.0 / (m_par.α - 1.0))
Kmin      = 1.0 * ((m_par.δ_0 - 0.0005 + (1.0 - m_par.β) / m_par.β) / m_par.α)^(0.5 / (m_par.α - 1.0))
K         = range(Kmin, stop = Kmax, length = 8)
# a.) Define excess demand function
d(K)      = Kdiff(K,n_par,m_par)

# b.) Find equilibrium capital stock (multigrid on y)
# ba.) initial calculation
KSS       = Brent(d, Kmin, Kmax)[1]

# c.) Calculate other equilibrium quantities
NSS       = employment(KSS, 1.0 ./ (m_par.μ*m_par.μw), m_par)
rSS       = interest(KSS,1.0 / m_par.μ, NSS, m_par)
wSS       = wage(KSS,1.0 / m_par.μ, NSS , m_par)
YSS       = output(KSS,1.0,NSS, m_par)
ProfitsSS = (1.0 -1.0 / m_par.μ).*YSS

KSS, BSS, TransitionMatSS, TransitionMat_aSS, TransitionMat_nSS, distrSS,
        c_a_starSS, m_a_starSS, k_a_starSS, c_n_starSS, m_n_starSS,VmSS, VkSS =
        Ksupply(m_par.RB./m_par.π,1.0+ rSS,wSS*NSS/n_par.H,ProfitsSS,n_par,m_par)

## bb.) refinement
ny              = n_par.ny_refined; # eigs in Ksupply quickly increases in runtime in ny
grid_y, Π, bounds= Tauchen(m_par.ρ_h,ny)
gy              = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2)); (m_par.ζ .+
                    m_par.ι)/m_par.ζ]
# Interpolate distribution on refined wealth-income grid
refined_dist    = mylinearinterpolate3(n_par.grid_m,n_par.grid_k, n_par.grid_y,
                                        distrSS,n_par.grid_m,n_par.grid_k, gy)
refined_dist    = refined_dist ./ sum(refined_dist,dims=(1,2,3))
# Include entrepreneurs into the income transitions
Π               = [Π .* (1.0 .- m_par.ζ)  m_par.ζ .* ones(ny);
                    m_par.ι ./ ny * ones(1,ny) 1.0 .- m_par.ι]
grid_y          = [exp.(grid_y .* m_par.σ_h ./ sqrt(1.0 .- m_par.ρ_h.^2));
                   (m_par.ζ .+ m_par.ι)/m_par.ζ]
Paux            = Π^1000
H               = (Paux[1,1:end-1]'*grid_y[1:end-1])

# Write changed parameter values to n_par
@set! n_par.ny          = ny + 1
@set! n_par.nstates     = (ny + 1) + n_par.nk + n_par.nm + length(state_names) - 3
@set! n_par.naggrstates = length(state_names)
@set! n_par.naggrcontrols = length(control_names)
@set! n_par.aggr_names  = aggr_names
@set! n_par.naggr       = length(aggr_names)

@set! n_par.bounds_y    = bounds
@set! n_par.ϵ           = 1e-10
@set! n_par.grid_y      = grid_y
@set! n_par.mesh_y      = repeat(reshape(grid_y, (1, 1, ny + 1)), outer=[n_par.nm, n_par.nk, 1])
@set! n_par.mesh_m      = repeat(reshape(n_par.grid_m, (n_par.nm, 1, 1)),
                        outer=[1, n_par.nk, ny + 1])
@set! n_par.mesh_k      = repeat(reshape(n_par.grid_k, (1, n_par.nk, 1)),
                                outer=[n_par.nm, 1, ny + 1])
@set! n_par.Π           = Π
@set! n_par.H           = H
@set! n_par.dist_guess  = refined_dist

KSS                     = Brent(d, KSS*.9, KSS*1.2)[1]

# c.) Calculate other equilibrium quantities
NSS                     = employment(KSS, 1.0 ./ (m_par.μ*m_par.μw), m_par)
rSS                     = interest(KSS,1.0 / m_par.μ, NSS, m_par)
wSS                     = wage(KSS,1.0 / m_par.μ, NSS , m_par)
YSS                     = output(KSS,1.0,NSS, m_par)
ProfitsSS               = (1.0 -1.0 / m_par.μ).*YSS

KSS, BSS, TransitionMatSS, TransitionMat_aSS, TransitionMat_nSS, distrSS,
        c_a_starSS, m_a_starSS, k_a_starSS, c_n_starSS, m_n_starSS,VmSS, VkSS =
        Ksupply(m_par.RB./m_par.π,1.0+ rSS,wSS*NSS/n_par.H,ProfitsSS,n_par,m_par)

VmSS                    = log.(VmSS)
VkSS                    = log.(VkSS)
RBSS                    = m_par.RB
ISS                     = m_par.δ_0*KSS

# Produce distributional summary statistics
eff_int      = (RBSS           .+ (m_par.Rbar .* (n_par.mesh_m.<=0.0))) # effective rate (need to check timing below and inflation)

incgross =[  (n_par.mesh_y .* (1. /m_par.μw).*wSS.*NSS./H).+
        (1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS,# labor income (NEW)
        interest(KSS,1.0 / m_par.μ, NSS, m_par).* n_par.mesh_k, # rental income
        eff_int .* n_par.mesh_m, # liquid asset Income
        n_par.mesh_k,
        (1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./H] # capital liquidation Income (q=1 in steady state)
incgross[1][:,:,end].= n_par.mesh_y[:,:,end] .* ProfitsSS  # profit income net of taxes
incgross[5][:,:,end].= n_par.mesh_y[:,:,end] .* ProfitsSS  # profit income net of taxes

inc =[  (((m_par.γ - m_par.τ_prog)/(m_par.γ+1)).*m_par.τ_lev.*(n_par.mesh_y.*1.0 ./m_par.μw.*wSS.*NSS./H).^(1.0-m_par.τ_prog)).+
        ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS),# labor income (NEW)
        interest(KSS,1.0 / m_par.μ, NSS, m_par).* n_par.mesh_k, # rental income
        eff_int .* n_par.mesh_m, # liquid asset Income
        n_par.mesh_k,
        m_par.τ_lev.*((1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./H).^(1.0-m_par.τ_prog).*((1.0 - m_par.τ_prog)/(m_par.γ+1)),
        m_par.τ_lev.*((1.0 ./ m_par.μw).*wSS.*NSS.*n_par.mesh_y./H).^(1.0-m_par.τ_prog)] # capital liquidation Income (q=1 in steady state)
inc[1][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)  # profit income net of taxes
inc[5][:,:,end].= 0.0
inc[6][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)  # profit income net of taxes

taxrev        = incgross[5]-inc[6]
incgrossaux   = incgross[5]
av_tax_rateSS = (distrSS[:]' * taxrev[:])./(distrSS[:]' * incgrossaux[:])

# apply taxes to union profits
inc[1] =(((m_par.γ - m_par.τ_prog)/(m_par.γ+1)).*m_par.τ_lev.*(n_par.mesh_y.*1.0 ./m_par.μw.*wSS.*NSS./H).^(1.0-m_par.τ_prog)).+ # labor income
        ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS).*(1.0 .- av_tax_rateSS)# labor union income
inc[1][:,:,end].= m_par.τ_lev.*(n_par.mesh_y[:,:,end] .* ProfitsSS).^(1.0-m_par.τ_prog)


TSS           = (distrSS[:]' * taxrev[:] + av_tax_rateSS*((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS))
GSS           = TSS - (m_par.RB./m_par.π-1.0)*BSS

distr_m_SS, distr_k_SS, distr_y_SS, share_borrowerSS, GiniWSS, I90shareSS,I90sharenetSS, GiniXSS,
        sdlogxSS, P9010CSS, GiniCSS, sdlgCSS, P9010ISS, GiniISS, sdlgISS, w90shareSS, P10CSS, P50CSS, P90CSS =
        distrSummaries(distrSS, c_a_starSS, c_n_starSS, n_par, inc,incgross, m_par)
# ------------------------------------------------------------------------------
## STEP 2: Dimensionality reduction
# ------------------------------------------------------------------------------
# 2a.) Discrete cosine transformation of policies or MUs
aux             = dct(VmSS)
ThetaVm   = aux[:] # Discrete Cosine transformation
ind             = sortperm(abs.(ThetaVm[:]);rev=true) #Indexes of sorted coefficients
coeffs          = 1
# Find the important basis functions (discrete cosine) for c_polSS
while norm(ThetaVm[ind[1:coeffs]])/norm(ThetaVm ) < 1 - n_par.reduc
        coeffs += 1
end
compressionIndexesVm  = ind[1:coeffs]

ThetaVk   = dct(VkSS)[:] # Discrete Cosine transformation
ind             = sortperm(abs.(ThetaVk[:]);rev=true) #Indexes of sorted coefficients
coeffs          = 1;
# Find the important basis functions (discrete cosine) for c_polSS
while norm(ThetaVk[ind[1:coeffs]])/norm(ThetaVk ) < 1 - n_par.reduc
        coeffs += 1
end
compressionIndexesVk  = ind[1:coeffs]

compressionIndexes =Array{Array{Int,1},1}(undef ,2)
compressionIndexes[1] = compressionIndexesVm
compressionIndexes[2] = compressionIndexesVk
# 2b.) Produce the Copula as an interpolant on the distribution function
#      and its marginals
CDF_SS     = zeros(n_par.nm+1,n_par.nk+1,n_par.ny+1)
CDF_SS[2:end,2:end,2:end]     = cumsum(cumsum(cumsum(distrSS,dims=1),dims=2),dims=3)
distr_m_SS = sum(distrSS,dims=(2,3))[:]
distr_k_SS = sum(distrSS,dims=(1,3))[:]
distr_y_SS = sum(distrSS,dims=(1,2))[:]
CDF_m      = cumsum([0.0; distr_m_SS[:]])
CDF_k      = cumsum([0.0; distr_k_SS[:]])
CDF_y      = cumsum([0.0; distr_y_SS[:]])

Copula(x::Vector,y::Vector,z::Vector) = mylinearinterpolate3(CDF_m, CDF_k, CDF_y,
                                                             CDF_SS, x, y, z)

# ------------------------------------------------------------------------------

@include "../input_aggregate_steady_state.jl"

# write to XSS vector
@writeXSS

# produce indexes to access XSS etc.
indexes = produce_indexes(n_par, compressionIndexesVm, compressionIndexesVk)
indexes_aggr = produce_indexes_aggr(n_par)

return XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, #=
        =# m_par, d, CDF_SS, CDF_m, CDF_k, CDF_y, distrSS
end
