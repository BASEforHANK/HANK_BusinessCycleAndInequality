
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
    # ------------------------------------------------------------------------------
    # STEP 1: Evaluate StE to calculate steady state variable values
    # ------------------------------------------------------------------------------
    # Calculate other equilibrium quantities
    Paux      = n_par.Π^1000                                              # Calculate ergodic ince distribution from transitions
    distr_y   = Paux[1, :] 
    NSS       = employment(KSS, 1.0./(m_par.μ*m_par.μw), m_par)
    rSS       = interest(KSS, 1.0./m_par.μ, NSS, m_par) + 1.0 
    wSS       = wage(KSS, 1.0./m_par.μ, NSS, m_par)
    YSS       = output(KSS, 1.0, NSS, m_par)                                               # stationary income distribution
    
    profitsSS       = (1.0 .- 1.0./m_par.μ) .*YSS
    unionprofitsSS  = (1.0 .- 1.0/m_par.μw) .* wSS .* NSS
    LC              = 1.0./m_par.μw *wSS.*NSS  
    taxrev          = ((n_par.grid_y/n_par.H).*LC)-m_par.τ_lev.*((n_par.grid_y/n_par.H).*LC).^(1.0-m_par.τ_prog)
    taxrev[end]     =  n_par.grid_y[end].*profitsSS - m_par.τ_lev.*( n_par.grid_y[end].*profitsSS).^(1.0-m_par.τ_prog)
    incgrossaux     = ((n_par.grid_y/n_par.H).*LC)
    incgrossaux[end]=  n_par.grid_y[end].*profitsSS
    av_tax_rateSS     = dot(distr_y, taxrev)./(dot(distr_y,incgrossaux))
    println(av_tax_rateSS)
    incgross, incnet, eff_int = 
            incomes(n_par, m_par, 1.0 ./ m_par.μw, 1.0, 1.0, 
                    m_par.RB, m_par.τ_prog, m_par.τ_lev, n_par.H, 1.0, 1.0,rSS,wSS,NSS,profitsSS,unionprofitsSS, av_tax_rateSS)
    # obtain other steady state variables
    KSS, BSS, TransitionMatSS, TransitionMat_aSS, TransitionMat_nSS,
                c_a_starSS, m_a_starSS, k_a_starSS, c_n_starSS, m_n_starSS, VmSS, VkSS, distrSS =
                        Ksupply(m_par.RB, rSS, n_par, m_par, VmSS, VkSS, distrSS, incnet, eff_int)
    
    VmSS        = log.(invmutil.(VmSS))
    VkSS        = log.(invmutil.(VkSS))
    
    
    # Produce distributional summary statistics
    distr_m_SS, distr_k_SS, distr_y_SS, TOP10WshareSS, TOP10IshareSS,TOP10InetshareSS, GiniWSS, GiniXSS, GiniCSS, GiniISS, sdlogySS = distrSummaries(distrSS, 1.0, c_a_starSS, c_n_starSS, n_par, incnet,incgross, m_par)
    
    # ------------------------------------------------------------------------------
    # STEP 2: Dimensionality reduction
    # ------------------------------------------------------------------------------
    # 2 a.) Discrete cosine transformation of marginal value functions
    # ------------------------------------------------------------------------------
    compressionIndexesVm  = select_comp_ind(VmSS,n_par.reduc_value) # store indexes of retained coefficients for Vm
    compressionIndexesVk  = select_comp_ind(VkSS,n_par.reduc_value) # store indexes of retained coefficients for Vk
    # ------------------------------------------------------------------------------
    # 2b.) Select polynomials for copula perturbation
    # ------------------------------------------------------------------------------
    SELECT                = [ ((i+j+k)<=n_par.reduc_copula) & 
                             (!((i==1) & (j ==1)) & !((k==1) & (j ==1)) & !((k==1) & (i ==1))) 
                             for i = 1:n_par.nm_copula, j = 1:n_par.nk_copula, k = 1:n_par.ny_copula]
    
    compressionIndexesCOP   = findall(SELECT[:])     # store indices of selected coeffs 

    # ------------------------------------------------------------------------------
    # 2c.) Store Compression Indexes
    # ------------------------------------------------------------------------------
    compressionIndexes    = Array{Array{Int,1},1}(undef ,3)       # Container to store all retained coefficients in one array
    compressionIndexes[1] = compressionIndexesVm
    compressionIndexes[2] = compressionIndexesVk
    compressionIndexes[3] = compressionIndexesCOP
    
    
    # ------------------------------------------------------------------------------
    # 2d.) Produce marginals
    # ------------------------------------------------------------------------------
    CDF_SS               = cumsum(cumsum(cumsum(distrSS,dims=1),dims=2),dims=3) # Calculate CDF from PDF
    distr_m_SS           = sum(distrSS,dims=(2,3))[:]            # Marginal distribution (pdf) of liquid assets
    distr_k_SS           = sum(distrSS,dims=(1,3))[:]            # Marginal distribution (pdf) of illiquid assets
    distr_y_SS           = sum(distrSS,dims=(1,2))[:]            # Marginal distribution (pdf) of income
    CDF_m                = cumsum(distr_m_SS[:])          # Marginal distribution (cdf) of liquid assets
    CDF_k                = cumsum(distr_k_SS[:])          # Marginal distribution (cdf) of illiquid assets
    CDF_y                = cumsum(distr_y_SS[:])          # Marginal distribution (cdf) of income
  
    function copula_marg_equi_y(distr_i,grid_i, nx)
       
        CDF_i                = cumsum(distr_i[:])          # Marginal distribution (cdf) of liquid assets

        distr_i_aux = distr_i[1:end-1]

        grid_i_aux = grid_i[1:end-1]

        aux_marginal      = log.(log.(collect(range(exp(exp(CDF_i[1])), stop=exp(exp(CDF_i[end])), length = nx))))

        function equishares(x1,x2) 
            
            FN_Wshares = cumsum(grid_i_aux.*distr_i_aux)./sum(grid_i_aux.*distr_i_aux);
            Wshares        = diff( HANKEstim.mylinearinterpolate(cumsum(distr_i_aux),FN_Wshares,[x1;x2]))
            loss=(Wshares .- 1.0/(nx-1))

        return loss
        end

        x2=1.0
        for i=2:nx-1
            equi(x1) = equishares(x1,x2) 

            x2 = find_zero(equi,(1e-9,x2))

            aux_marginal[end-i] = x2

        end

        aux_marginal[end-1] = CDF_i[end-1]
        aux_marginal[end] = CDF_i[end]
        aux_marginal[1]   = CDF_i[1]
        copula_marginal   = copy(aux_marginal)
        jlast             = nx
        for i = nx-1:-1:1
            j = HANKEstim.locate(aux_marginal[i], CDF_i) + 1
            if jlast == j 
                j -=1
            end
            jlast = j
            copula_marginal[i] = CDF_i[j]
        end
        return copula_marginal
    end

    function copula_marg_equi(distr_i,grid_i, nx)

        CDF_i                = cumsum(distr_i[:])          # Marginal distribution (cdf) of liquid assets

        aux_marginal      = log.(log.(collect(range(exp(exp(CDF_i[1])), stop=exp(exp(CDF_i[end])), length = nx))))

        function equishares(x1,x2) 
            
            FN_Wshares = cumsum(grid_i.*distr_i)./sum(grid_i.*distr_i);
            Wshares        = diff( HANKEstim.mylinearinterpolate(cumsum(distr_i),FN_Wshares,[x1;x2]))
            loss=(Wshares .- 1.0/(nx-1))

        return loss
        end

        x2=1.0
        for i=1:nx-1
            equi(x1) = equishares(x1,x2) 

            x2 = find_zero(equi,(1e-9,x2))

            aux_marginal[end-i] = x2

        end

        aux_marginal[end] = CDF_i[end]
        aux_marginal[1]   = CDF_i[1]
        copula_marginal   = copy(aux_marginal)
        jlast             = nx
        for i = nx-1:-1:1
            j = HANKEstim.locate(aux_marginal[i], CDF_i) + 1
            if jlast == j 
                j -=1
            end
            jlast = j
            copula_marginal[i] = CDF_i[j]
        end
        return copula_marginal
    end

    @set! n_par.copula_marginal_m = copula_marg_equi(distr_m_SS,n_par.grid_m,n_par.nm_copula)
    @set! n_par.copula_marginal_k = copula_marg_equi(distr_k_SS,n_par.grid_k,n_par.nk_copula)
    @set! n_par.copula_marginal_y = copula_marg_equi_y(distr_y_SS,n_par.grid_y,n_par.ny_copula)

    # ------------------------------------------------------------------------------
    
    # aggregate steady state marker
    # @include "../3_Model/input_aggregate_steady_state.jl"
    
    # write to XSS vector
    @writeXSS
    
    # produce indexes to access XSS etc.
    indexes               = produce_indexes(n_par, compressionIndexesVm, compressionIndexesVk, compressionIndexesCOP)
    indexes_aggr          = produce_indexes_aggr(n_par)
   
    @set! n_par.ntotal    = length(vcat(compressionIndexes...)) + (n_par.ny + n_par.nm + n_par.nk - 3 + n_par.naggr) 
    @set! n_par.nstates   = n_par.ny + n_par.nk + n_par.nm - 3 + n_par.naggrstates + length(compressionIndexes[3]) # add to no. of states the coefficients that perturb the copula
    @set! n_par.ncontrols = length(vcat(compressionIndexes[1:2]...)) + n_par.naggrcontrols
    @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
    @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
    
    if n_par.n_agg_eqn != n_par.naggr - length(n_par.distr_names)
        @warn("Inconsistency in number of aggregate variables and equations")
    end

    return XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, n_par, #=
            =# m_par, CDF_SS, CDF_m, CDF_k, CDF_y, distrSS       
    end