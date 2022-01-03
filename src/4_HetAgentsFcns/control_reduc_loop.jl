function control_reduc_loop(compressionIndexes, XSS, m_par, n_par, indexes, distrSS, shock_names)

  #---------------------------------------------------------------
  ## Step 1: Solve model around parameter guess
  #---------------------------------------------------------------
  A = zeros(n_par.ntotal, n_par.ntotal)
  B = zeros(n_par.ntotal, n_par.ntotal)
  State2Control, LOMstate, SolutionError, nk = SGU(XSS, A, B, m_par, n_par, indexes, compressionIndexes, distrSS; estim = false);  
 
  #---------------------------------------------------------------
  ## STEP 2: PRODUCE IRFs
  #---------------------------------------------------------------
  nlag = 200
  n_shocks = length(shock_names)
  IRF_state_sparse = zeros(indexes.profits, nlag, n_shocks)
  j = 0
  for i in shock_names
      x0      = zeros(size(LOMstate, 1), 1)
      j += 1
      x0[getfield(indexes, i)]  = getfield(m_par, Symbol("σ_", i))
      MX = [I; State2Control]
      x_aux = copy(x0)
      for t = 1:nlag
          IRF_state_sparse[:, t, j] = (MX*x_aux)'
          x_aux[:] = LOMstate * x_aux
      end
  end
  VARdecomp = copy(IRF_state_sparse)
  for i = 1:n_shocks
      VARdecomp[:, :, i] = cumsum(IRF_state_sparse[:, :, i].^2.0, dims=2)
  end
  
  #-------------------------------------------------------------
  ## Step 3: loop over horizons and retain only most volatile coefficients
  #-------------------------------------------------------------

  irf_horizon = 200
  reduc_crit  = 0.999
  
  IXVmtotal   = red_loop(VARdecomp[indexes.Vm,:,:], irf_horizon,reduc_crit)
  compressionIndexesVm_red  = compressionIndexes[1][IXVmtotal]
  compressionIndexes[1]     = compressionIndexesVm_red

  IXVktotal   = red_loop(VARdecomp[indexes.Vk,:,:], irf_horizon,reduc_crit)
  compressionIndexesVk_red  = compressionIndexes[2][IXVktotal] 
  compressionIndexes[2]     = compressionIndexesVk_red

#   reduc_crit  = 0.9
#   IXCOPtotal  = red_loop(VARdecomp[indexes.COP,:,:], irf_horizon,reduc_crit)
#   compressionIndexesCOP_red = compressionIndexes[3][IXCOPtotal]
#   compressionIndexes[3]     = compressionIndexesCOP_red

  indexes                        = produce_indexes(n_par, compressionIndexes[1], compressionIndexes[2], compressionIndexes[3])

  @set! n_par.ntotal    = length(vcat(compressionIndexes...)) + (n_par.ny + n_par.nm + n_par.nk - 3 + n_par.naggr) 
  @set! n_par.nstates   = n_par.ny + n_par.nk + n_par.nm - 3 + n_par.naggrstates + length(compressionIndexes[3]) # add to no. of states the coefficients that perturb the copula
  @set! n_par.ncontrols = length(vcat(compressionIndexes[1:2]...)) + n_par.naggrcontrols
  @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
  @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
   
  return compressionIndexes, indexes, n_par
end

function red_loop(VCD, irf_horizon,reduc_crit)
    AA      = sum(VCD[:,1,:], dims = 2)
    BB      = AA./sum(AA)
    IX      = sortperm(BB[:], rev=true)
    VAR     = cumsum(BB[IX])
    keep    = max(count(VAR .< reduc_crit), 1) 
    IXtotal = IX[1:keep]

    for tt = 2:irf_horizon
        AA      = sum(VCD[:,tt,:], dims = 2)
        BB      = AA./sum(AA)
        IX      = sortperm(BB[:], rev=true)
        count   = 0
        energy  = sum(BB[IXtotal])
        while energy < reduc_crit
              count  += 1
              IXtotal = union(IXtotal,IX[count])
              energy  = sum(BB[IXtotal])
        end     
    end
    return IXtotal
end