function MakeTransition(m_a_star::Array{Float64,3},
    m_n_star::Array{Float64,3},
    k_a_star::Array{Float64,3},
    Π::Array{Float64,2}, n_par::NumericalParameters)

    # create linear interpolation weights from policy functions
    idk_a, weightright_k_a, weightleft_k_a = MakeWeights(k_a_star,n_par.grid_k)
    idm_a, weightright_m_a, weightleft_m_a = MakeWeights(m_a_star,n_par.grid_m)
    idm_n, weightright_m_n, weightleft_m_n = MakeWeights(m_n_star,n_par.grid_m)

    # Adjustment case
    weight      = Array{typeof(k_a_star[1]),3}(undef, 4,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    targetindex = zeros(Int,4,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    startindex  = zeros(Int,4,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    blockindex  = (0:n_par.ny-1)*n_par.nk*n_par.nm
    runindex    = 0

    for zz = 1:n_par.ny # all current income states
        for kk = 1:n_par.nk # all current illiquid asset states
            for mm = 1:n_par.nm # all current liquid asset states
                runindex = runindex+1
                WLL      = weightleft_m_a[mm,kk,zz] .* weightleft_k_a[mm,kk,zz]
                WRL      = weightright_m_a[mm,kk,zz].* weightleft_k_a[mm,kk,zz]
                WLR      = weightleft_m_a[mm,kk,zz] .* weightright_k_a[mm,kk,zz]
                WRR      = weightright_m_a[mm,kk,zz].* weightright_k_a[mm,kk,zz]
                IDD      = idm_a[mm,kk,zz].+(idk_a[mm,kk,zz]-1).*n_par.nm
                for jj = 1:n_par.ny
                    pp                         = Π[zz,jj]
                    bb                         = blockindex[jj]
                    weight[1,jj,runindex]      = WLL .* pp
                    weight[2,jj,runindex]      = WRL .* pp
                    weight[3,jj,runindex]      = WLR .* pp
                    weight[4,jj,runindex]      = WRR .* pp
                    targetindex[1,jj,runindex] = IDD .+ bb
                    targetindex[2,jj,runindex] = IDD + 1 .+ bb
                    targetindex[3,jj,runindex] = IDD + n_par.nm .+ bb
                    targetindex[4,jj,runindex] = IDD + n_par.nm + 1 .+ bb
                    startindex[1,jj,runindex]  = runindex
                    startindex[2,jj,runindex]  = runindex
                    startindex[3,jj,runindex]  = runindex
                    startindex[4,jj,runindex]  = runindex
                end
            end
        end
    end
    S_a          = startindex[:]
    T_a          = targetindex[:]
    W_a          = weight[:]

    # Non-Adjustment case
    weight2      = zeros(typeof(k_a_star[1]), 2,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    targetindex2 = zeros(Int, 2,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    startindex2  = zeros(Int,2,n_par.ny,n_par.nk* n_par.nm*n_par.ny)
    runindex     = 0
    for zz = 1:n_par.ny # all current income states
        for kk = 1:n_par.nk # all current illiquid asset states
            for mm = 1:n_par.nm # all current liquid asset states
                runindex = runindex+1
                WL       = weightleft_m_n[mm,kk,zz]
                WR       = weightright_m_n[mm,kk,zz]
                CI       = idm_n[mm,kk,zz].+(kk-1).*n_par.nm
                for jj = 1:n_par.ny
                    pp                          = Π[zz,jj]
                    weight2[1,jj,runindex]      = WL .* pp
                    weight2[2,jj,runindex]      = WR .* pp
                    targetindex2[1,jj,runindex] = CI .+ blockindex[jj]
                    targetindex2[2,jj,runindex] = CI .+ 1 .+blockindex[jj]
                    startindex2[1,jj,runindex]  = runindex
                    startindex2[2,jj,runindex]  = runindex
                end
            end
        end
    end
    S_n        = startindex2[:]
    T_n        = targetindex2[:]
    W_n        = weight2[:]

    return S_a, T_a, W_a, S_n, T_n, W_n
end
