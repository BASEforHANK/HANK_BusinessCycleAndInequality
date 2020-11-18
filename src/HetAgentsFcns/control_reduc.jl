function control_reduc(compressionIndexes, XSS, m_par, n_par, indexes, Copula, distrSS, shock_names)

  #---------------------------------------------------------------
  ## Step 1: Solve model around parameter guess
  #---------------------------------------------------------------
  A = zeros(n_par.ntotal,n_par.ntotal)
  B = zeros(n_par.ntotal,n_par.ntotal)
  State2Control, LOMstate, SolutionError, nk = SGU(XSS, A, B, m_par, n_par, indexes, Copula, compressionIndexes, distrSS; estim = false);  
 
  #---------------------------------------------------------------
  ## STEP 2: PRODUCE IRFs
  #---------------------------------------------------------------
  nlag = 1000
  n_shocks = length(shock_names)
  IRF_state_sparse = zeros(indexes.profits, nlag, n_shocks)
  j=0
  for i in shock_names
      x0      = zeros(size(LOMstate,1), 1)
      j += 1
      x0[getfield(indexes,i)]  = getfield(m_par,Symbol("σ_",i))
      MX = [I; State2Control]
      x_aux = copy(x0)
      for t = 1:nlag
          IRF_state_sparse[:, t, j]= (MX*x_aux)'
          x_aux[:] = LOMstate * x_aux
      end
  end
  VARdecomp = copy(IRF_state_sparse)
  for i=1:n_shocks
      VARdecomp[:,:,i] = cumsum(IRF_state_sparse[:,:,i].^2.0, dims=2)
  end
  #-------------------------------------------------------------
  ## Step 3: Variance decomposition at long horizon
  #-------------------------------------------------------------
  CVD = VARdecomp[:,end,:];
  AA  = sum(CVD[indexes.Vk,:],dims=2)
  BB1 = AA./sum(AA)
  AA  = sum(CVD[indexes.Vm,:],dims=2)
  BB2 = AA./sum(AA)
  #-------------------------------------------------------------
  ## Step 4: Select basis functions important at long horizon
  #-------------------------------------------------------------
  IXVk   = sortperm(BB1[:],rev=true )
  VARVk  = cumsum(BB1[IXVk])
  keepVk = count(VARVk.<0.9999)+1
  IXVktotal = IXVk[1:keepVk]

  IXVm   = sortperm(BB2[:],rev=true)
  VARVm  = cumsum(BB2[IXVm])
  keepVm = count(VARVm.<0.9999)+1
  IXVmtotal=IXVm[1:keepVm]

  #-------------------------------------------------------------
  ## Step 5: Add basis functions of particular importance at 
  ##         short horizons
  #------------------------------------------------------------- 

  for tt=1:1000-1
      CVD=VARdecomp[:,tt,:];
      AA=sum(CVD[indexes.Vk,:],dims=2)
      BB1=AA./sum(AA)

      AA=sum(CVD[indexes.Vm,:],dims=2)
      BB2=AA./sum(AA)

      IXVk2=sortperm(BB1[:],rev=true)
      VARVk2=cumsum(BB1[IXVk2])
      keepVk2=count(VARVk2.<0.999)+1
      IXVktotal=unique(vcat(IXVktotal,IXVk2[1:keepVk2]))

      IXVm2=sortperm(BB2[:],rev=true)
      VARVm2 = cumsum(BB2[IXVm2])
      keepVm2=count(VARVm2.<0.999)+1

      IXVmtotal=unique(vcat(IXVmtotal,IXVm2[1:keepVm2] ))
  end

  keepVmtotal = length(IXVmtotal)

  keepVktotal = length(IXVktotal)

  compressionIndexesVm_red = compressionIndexes[1][IXVmtotal]
  compressionIndexesVk_red = compressionIndexes[2][IXVktotal]
  compressionIndexes[1] = compressionIndexesVm_red
  compressionIndexes[2] = compressionIndexesVk_red

  indexes      = produce_indexes(n_par, compressionIndexes[1], compressionIndexes[2], compressionIndexes[3])
  ntotal                  = indexes.profits    
  @set! n_par.ntotal      = ntotal
  @set! n_par.ncontrols   = length(compressionIndexes[1]) + length(compressionIndexes[2]) + n_par.naggrcontrols
  @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
  @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)

  return compressionIndexes, indexes, n_par
end
