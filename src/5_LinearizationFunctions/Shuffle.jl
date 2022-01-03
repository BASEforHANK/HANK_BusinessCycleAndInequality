##########################################################
# Matrix to remove one degree of freedom from distribution
#---------------------------------------------------------
function shuffleMatrix(distr, n_par)
    distr_m = sum(sum(distr,dims=3),dims=2)./sum(distr[:])
    distr_k = sum(sum(distr,dims=3),dims=1)./sum(distr[:])
    distr_y = sum(sum(distr,dims=2),dims=1)./sum(distr[:])
    Γ       = Array{Array{Float64,2},1}(undef,3)
    Γ[1]    = zeros(Float64,(n_par.nm,n_par.nm-1))
    Γ[2]    = zeros(Float64,(n_par.nk,n_par.nk-1))
    Γ[3]    = zeros(Float64,(n_par.ny,n_par.ny-1))
    for j=1:n_par.nm-1
        Γ[1][:,j] = -distr_m[:]
        Γ[1][j,j] = 1-distr_m[j]
        Γ[1][j,j] = Γ[1][j,j] - sum(Γ[1][:,j])
    end
    for j=1:n_par.nk-1
        Γ[2][:,j] = -distr_k[:]
        Γ[2][j,j] = 1-distr_k[j]
        Γ[2][j,j] = Γ[2][j,j] - sum(Γ[2][:,j])
    end
    for j=1:n_par.ny-1
        Γ[3][:,j] = -distr_y[:]
        Γ[3][j,j] = 1-distr_y[j]
        Γ[3][j,j] = Γ[3][j,j] - sum(Γ[3][:,j])
    end

    return Γ
end
