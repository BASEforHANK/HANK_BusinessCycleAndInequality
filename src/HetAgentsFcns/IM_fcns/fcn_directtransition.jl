function DirectTransition(m_a_star::Array,
    m_n_star::Array,
    k_a_star::Array,
    distr::Array,
    λ,
    Π::Array,
    n_par::NumericalParameters)

    dPrime = zeros(eltype(distr),size(distr))
    idk_a, wR_k_a = MakeWeightsLight(k_a_star,n_par.grid_k)
    idm_a, wR_m_a = MakeWeightsLight(m_a_star,n_par.grid_m)
    idm_n, wR_m_n = MakeWeightsLight(m_n_star,n_par.grid_m)
    blockindex = (0:n_par.ny-1)*n_par.nk*n_par.nm
    @views @inbounds begin
    for zz = 1:n_par.ny # all current income states
        for kk = 1:n_par.nk # all current illiquid asset states
            #idk_n = kk
            for mm = 1:n_par.nm
                dd=distr[mm,kk,zz]
                IDD_a = (idm_a[mm,kk,zz].+(idk_a[mm,kk,zz] .-1).*n_par.nm)
                IDD_n = (idm_n[mm,kk,zz].+(kk-1).*n_par.nm)
                dl    = (dd.*(1.0 .- wR_k_a[mm,kk,zz]))
                DLL_a = (dl.*(1.0 .- wR_m_a[mm,kk,zz]))
                DLR_a = (dl.*wR_m_a[mm,kk,zz])
                dr    = (dd.*wR_k_a[mm,kk,zz])
                DRL_a = (dr.*(1.0 .- wR_m_a[mm,kk,zz]))
                DRR_a = (dr.*wR_m_a[mm,kk,zz])
                DL_n  = (dd.*(1.0 .- wR_m_n[mm,kk,zz]))
                DR_n  = (dd.*wR_m_n[mm,kk,zz])
                pp    = (Π[zz,:])
                for yy = 1:n_par.ny
                    id_a = IDD_a .+ blockindex[yy]
                    id_n = IDD_n .+ blockindex[yy]
                    fac = λ.*pp[yy]
                    dPrime[id_a]            += fac.*DLL_a
                    dPrime[id_a+1]          += fac.*DLR_a
                    dPrime[id_a+n_par.nm]   += fac.*DRL_a
                    dPrime[id_a+n_par.nm+1] += fac.*DRR_a
                    dPrime[id_n]            += (1.0 .- λ).*pp[yy].*DL_n
                    dPrime[id_n+1]          += (1.0 .- λ).*pp[yy].*DR_n
                end
            end
        end
    end
end
    return dPrime
end
