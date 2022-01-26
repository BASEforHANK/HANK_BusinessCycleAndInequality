@doc raw"""
    compute_reduction(sr, lr, m_par, shock_names)

Finds the coefficients of the Chebyshev Polynomials obtained from the DCTs of marginal value functions and copulas that are most volatile over the business cycle.
Thereby, it allows to reduce the number of state and control variables required to describe the heterogeneous agent part of the model.
This improves the speed of estimation without losing much of the approximation quality.

# Returns
- `compressionIndexes`: retained DCT coefficients (stored in an array of arrays)
- `indexes`: indexes that translate variable names to positions in arrays
- `n_par`: numerical parameters (changed number of states and controls)
"""
function compute_reduction(sr, lr, m_par, shock_names)

  #---------------------------------------------------------------
  ## STEP 1: PRODUCE Long Run Covariance
  #---------------------------------------------------------------
  n_par = sr.n_par    
  SCov = zeros(n_par.nstates, n_par.nstates)
  for i in shock_names
      SCov[getfield(sr.indexes, i), getfield(sr.indexes, i)] = (getfield(m_par,Symbol("σ_", i))).^2
  end
  
  StateCOVAR   = lyapd(lr.LOMstate, SCov)
  
  ControlCOVAR = lr.State2Control*StateCOVAR*lr.State2Control'
  ControlCOVAR = (ControlCOVAR + ControlCOVAR')./2
  #---------------------------------------------------------------
  ## STEP 1: Produce eigenvalue decomposition of relevant COVARs
  #          and select eigenvectors of large eigenvalues to obtain 
  #          low-dim to high-dim projections of copula and value-functions
  #---------------------------------------------------------------

  Dindex       = sr.indexes.COP
  evalS, evecS = eigen(StateCOVAR[Dindex,Dindex])
  keepD        = abs.(evalS).> maximum(evalS)*n_par.further_compress_critS
  indKeepD     = Dindex[keepD]
  @set! n_par.nstates_r = n_par.nstates - length(Dindex) + length(indKeepD) 

  Vindex       = [sr.indexes.Vm; sr.indexes.Vk] 
  evalC, evecC = eigen(ControlCOVAR[Vindex.- n_par.nstates, Vindex.- n_par.nstates])
  keepV        = abs.(evalC).> maximum(evalC)*n_par.further_compress_critC
  indKeepV     = Vindex[keepV]
  @set! n_par.ncontrols_r = n_par.ncontrols - length(Vindex) + length(indKeepV) 
  
  #-------------------------------------------------------------
  ## Step 3: Put together projection matrices and update indexes
  #-------------------------------------------------------------

  PRightStates_aux = float(I[1:n_par.nstates, 1:n_par.nstates])
  PRightStates_aux[Dindex,Dindex] = evecS
  keep             = ones(Bool,n_par.nstates)
  keep[Dindex[.!keepD]]   .= false 
  @set! n_par.PRightStates = PRightStates_aux[:,keep]
    
  PRightAll_aux= float(I[1:n_par.ntotal, 1:n_par.ntotal])
  PRightAll_aux[Dindex,Dindex] = evecS
  PRightAll_aux[Vindex,Vindex] = evecC
  keep = ones(Bool,n_par.ntotal)
  keep[Dindex[.!keepD]]       .= false 
  keep[Vindex[.!keepV]]       .= false 
  @set! n_par.PRightAll        = PRightAll_aux[:,keep]

  indexes_r                    = produce_indexes(n_par, keepV[keepV][1:2], keepV[keepV][3:end], keepD[keepD]) # arbitrary splitup of underlying factors from value functions to indexes in reduced model

  @set! n_par.ntotal_r         = n_par.nstates_r + n_par.ncontrols_r 

  return  indexes_r, n_par
end