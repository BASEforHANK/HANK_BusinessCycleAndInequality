function MultipleDirectTransition(m_a_star::Array{Float64,3},
    m_n_star::Array{Float64,3},
    k_a_star::Array{Float64,3},
    distr::Array{Float64,3},
    λ::Float64,
    Π::Array{Float64,2}, n_par::NumericalParameters)


    idk_a, wR_k_a, wL_k_a = MakeWeights(k_a_star,n_par.grid_k)
    idm_a, wR_m_a, wL_m_a = MakeWeights(m_a_star,n_par.grid_m)
    idm_n, wR_m_n, wL_m_n = MakeWeights(m_n_star,n_par.grid_m)
    dist =1.0
    count = 1
    blockindex = (0:n_par.ny-1)*n_par.nk*n_par.nm
    while (dist>n_par.ϵ) && (count<10000)
        dPrime = zeros(typeof(distr[1]),size(distr))
        for zz = 1:n_par.ny # all current income states
            for kk = 1:n_par.nk # all current illiquid asset states
                #idk_n = kk
                for mm = 1:n_par.nm
                    dd=distr[mm,kk,zz]
                    IDD_a = idm_a[mm,kk,zz].+(idk_a[mm,kk,zz]-1).*n_par.nm
                    IDD_n = idm_n[mm,kk,zz].+(kk-1).*n_par.nm
                    DLL_a = dd.*wL_k_a[mm,kk,zz].*wL_m_a[mm,kk,zz]
                    DLR_a = dd.*wL_k_a[mm,kk,zz].*wR_m_a[mm,kk,zz]
                    DRL_a = dd.*wR_k_a[mm,kk,zz].*wL_m_a[mm,kk,zz]
                    DRR_a = dd.*wR_k_a[mm,kk,zz].*wR_m_a[mm,kk,zz]
                    DL_n  = dd.*wL_m_n[mm,kk,zz]
                    DR_n  = dd.*wR_m_n[mm,kk,zz]
                    for yy = 1:n_par.ny
                        id_a = IDD_a +blockindex[yy]
                        id_n = IDD_n +blockindex[yy]
                        fac = λ.*Π[zz,yy]
                        dPrime[id_a]            += fac.*DLL_a
                        dPrime[id_a+1]          += fac.*DLR_a
                        dPrime[id_a+n_par.nm]   += fac.*DRL_a
                        dPrime[id_a+n_par.nm+1] += fac.*DRR_a
                        dPrime[id_n]            += (1.0-λ).*Π[zz,yy].*DL_n
                        dPrime[id_n+1]          += (1.0-λ).*Π[zz,yy].*DR_n
                    end
                end
            end
        end
        dist  = maximum(abs.(dPrime[:] - distr[:]))
        distr = dPrime
        count = count+1
    end
    return distr, dist, count
end
