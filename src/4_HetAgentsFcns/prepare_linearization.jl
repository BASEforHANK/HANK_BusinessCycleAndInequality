
@doc raw"""
    prepare_linearization(KSS, VmSS, VkSS, distrSS, n_par, m_par)

Compute a number of equilibrium objects needed for linearization.

# Arguments
- `KSS`: steady-state capital stock
- `VmSS`, `VkSS`: marginal value functions
- `distrSS::Array{Float64,3}`: steady-state distribution of idiosyncratic states, computed by [`Ksupply()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`

# Returns
- `XSS::Array{Float64,1}`, `XSSaggr::Array{Float64,1}`: steady state vectors produced by [`@writeXSS()`](@ref)
- `indexes`, `indexes_aggr`: `struct`s for accessing `XSS`,`XSSaggr` by variable names, produced by [`@make_fn()`](@ref),
        [`@make_fnaggr()`](@ref)
- `compressionIndexes::Array{Array{Int,1},1}`: indexes for compressed marginal value functions (``V_m`` and ``V_k``)
- `Copula(x,y,z)`: function that maps marginals `x`,`y`,`z` to approximated joint distribution, produced by
        [`mylinearinterpolate3()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`
- `CDF_SS`, `CDF_m`, `CDF_k`, `CDF_y`: cumulative distribution functions (joint and marginals)
- `distrSS::Array{Float64,3}`: steady state distribution of idiosyncratic states, computed by [`Ksupply()`](@ref)
"""
function prepare_linearization(KSS, VmSS, VkSS, distrSS, n_par, m_par)
    if n_par.verbose
        println("Running reduction step to prepare linearization")
    end
    # Calculate other equilibrium quantities
    incgross, incnet, NSS, rSS, wSS, YSS, ProfitsSS, ISS, RBSS, taxrev, av_tax_rateSS, eff_int = incomes(n_par, m_par, KSS, distrSS)
    
    # obtain other steady state variables
    KSS, BSS, TransitionMatSS, TransitionMat_aSS, TransitionMat_nSS,
        c_a_starSS, m_a_starSS, k_a_starSS, c_n_starSS, m_n_starSS, VmSS, VkSS, distrSS =
                          Ksupply(RBSS, 1.0 + rSS, n_par, m_par, VmSS, VkSS, distrSS, incnet, eff_int)
    
    VmSS                = log.(VmSS)
    VkSS                = log.(VkSS)
    # Calculate taxes and government expenditures
    TSS                 = (distrSS[:]' * taxrev[:] + av_tax_rateSS*((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS))
    GSS                 = TSS - (m_par.RB./m_par.π-1.0)*BSS
    
    # Produce distributional summary statistics
    distr_m_SS, distr_k_SS, distr_y_SS, share_borrowerSS, GiniWSS, I90shareSS,I90sharenetSS, GiniXSS,
            sdlogxSS, P9010CSS, GiniCSS, sdlgCSS, P9010ISS, GiniISS, sdlogySS, w90shareSS, P10CSS, P50CSS, P90CSS =
            distrSummaries(distrSS, c_a_starSS, c_n_starSS, n_par, incnet,incgross, m_par)
    # ------------------------------------------------------------------------------
    ## STEP 2: Dimensionality reduction
    # ------------------------------------------------------------------------------
    # 2 a.) Discrete cosine transformation of marginal value functions
    # ------------------------------------------------------------------------------
    ThetaVm             = dct(VmSS)[:]                          # Discrete cosine transformation of marginal liquid asset value
    ind                 = sortperm(abs.(ThetaVm[:]);rev=true)   # Indexes of coefficients sorted by their absolute size
    coeffs              = 1                                     # Container to store the number of retained coefficients
    # Find the important basis functions (discrete cosine) for VmSS
    while norm(ThetaVm[ind[1:coeffs]])/norm(ThetaVm ) < 1 - n_par.reduc_value 
            coeffs      += 1                                    # add retained coefficients until only n_par.reduc_value hare of energy is lost
    end
    compressionIndexesVm  = ind[1:coeffs]                       # store indexes of retained coefficients
    
    ThetaVk             = dct(VkSS)[:]                          # Discrete cosine transformation of marginal illiquid asset value
    ind                 = sortperm(abs.(ThetaVk[:]);rev=true)   # Indexes of coefficients sorted by their absolute size
    coeffs              = 1                                     # Container to store the number of retained coefficients
    # Find the important basis functions (discrete cosine) for VkSS
    while norm(ThetaVk[ind[1:coeffs]])/norm(ThetaVk ) < 1 - n_par.reduc_value
            coeffs      += 1                                    # add retained coefficients until only n_par.reduc_value hare of energy is lost
    end
    compressionIndexesVk = ind[1:coeffs]                        # store indexes of retained coefficients
    
    distr_LOL           = distrSS[1:end-1,1:end-1,1:end-1]      # Leave out last entry of histogramm (b/c it integrates to 1)
    ThetaD              = dct(distr_LOL)[:]                     # Discrete cosine transformation of Copula
    ind                 = sortperm(abs.(ThetaD[:]);rev=true)    # Indexes of coefficients sorted by their absolute size
    compressionIndexesD = ind[2:2+n_par.reduc_copula]           # keep n_par.reduc_copula coefficients, leave out index no. 1 as this shifts the constant
    
    compressionIndexes  = Array{Array{Int,1},1}(undef ,3)       # Container to store all retained coefficients in one array
    compressionIndexes[1] = compressionIndexesVm
    compressionIndexes[2] = compressionIndexesVk
    compressionIndexes[3] = compressionIndexesD
    
    @set! n_par.nstates = n_par.ny + n_par.nk + n_par.nm + n_par.naggrstates - 3 + length(compressionIndexes[3]) # add to no. of states the coefficients that perturb the copula
    # ------------------------------------------------------------------------------
    # 2b.) Produce the Copula as an interpolant on the distribution function
    #      and its marginals
    # ------------------------------------------------------------------------------
    CDF_SS              = zeros(n_par.nm+1,n_par.nk+1,n_par.ny+1) # Produce CDF of asset-income distribution (container here)
    CDF_SS[2:end,2:end,2:end]     = cumsum(cumsum(cumsum(distrSS,dims=1),dims=2),dims=3) # Calculate CDF from PDF
    distr_m_SS          = sum(distrSS,dims=(2,3))[:]            # Marginal distribution (pdf) of liquid assets
    distr_k_SS          = sum(distrSS,dims=(1,3))[:]            # Marginal distribution (pdf) of illiquid assets
    distr_y_SS          = sum(distrSS,dims=(1,2))[:]            # Marginal distribution (pdf) of income
    CDF_m               = cumsum([0.0; distr_m_SS[:]])          # Marginal distribution (cdf) of liquid assets
    CDF_k               = cumsum([0.0; distr_k_SS[:]])          # Marginal distribution (cdf) of illiquid assets
    CDF_y               = cumsum([0.0; distr_y_SS[:]])          # Marginal distribution (cdf) of income
    
    Copula(x::Vector,y::Vector,z::Vector) = 
        mylinearinterpolate3(CDF_m, CDF_k, CDF_y,CDF_SS, x, y, z) # Define Copula as a function (used only if not perturbed)
    
    # ------------------------------------------------------------------------------
    
    # aggregate steady state marker
    # @include "../3_Model/input_aggregate_steady_state.jl"
    
    # write to XSS vector
    @writeXSS
    
    # produce indexes to access XSS etc.
    indexes               = produce_indexes(n_par, compressionIndexesVm, compressionIndexesVk, compressionIndexesD)
    indexes_aggr          = produce_indexes_aggr(n_par)
    ntotal                = indexes.profits # Convention: profits is the last control in the list of control variables
    @set! n_par.ntotal    = ntotal   
    @set! n_par.ncontrols = length(compressionIndexes[1]) + length(compressionIndexes[2]) + n_par.naggrcontrols
    @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
    @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
    
    return XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, #=
            =# m_par, CDF_SS, CDF_m, CDF_k, CDF_y, distrSS       
    end