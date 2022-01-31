function compare_2_linearizations(sr_reduc,lr_reduc,m_par_reduc,sr_full,lr_full, m_par_full, legend_entry)
    IRF_r = irfs(sr_reduc, lr_reduc,m_par_reduc)
    IRF_f = irfs(sr_full, lr_full,m_par_full)

    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 6)
    p[1] = plot(IRF_r[sr_reduc.indexes_r.Z, :], title = "TFP (percent)", label = false)
    p[1] = plot!(IRF_f[sr_full.indexes_r.Z, :], label = false, linestyle = :dash)
    p[2] = plot(IRF_r[sr_reduc.indexes_r.I, :], title = "Investment (percent)", label = false)
    p[2] = plot!(IRF_f[sr_full.indexes_r.I, :], label = false, linestyle = :dash)
    p[3] = plot(IRF_r[sr_reduc.indexes_r.Y, :], title = "GDP (percent)", label = false)
    p[3] = plot!(IRF_f[sr_full.indexes_r.Y, :], label = false, linestyle = :dash)
    p[4] = plot(IRF_r[sr_reduc.indexes_r.C, :], title = "Consumption (percent)", label = false)
    p[4] = plot!(IRF_f[sr_full.indexes_r.C, :], label = false, linestyle = :dash)
    p[5] = plot(IRF_r[sr_reduc.indexes_r.TOP10Wshare, :], title = "Top 10 wealth share (percent)", label = false)
    p[5] = plot!(IRF_f[sr_full.indexes_r.TOP10Wshare, :], linestyle = :dash, label = false)
    p[6] = plot(Matrix{Missing}(undef, 3,1),  label = legend_entry[1], framestyle=:none, legendfontsize=12)
    p[6] = plot!(Matrix{Missing}(undef, 3,1), label = legend_entry[2], framestyle=:none, linestyle = :dash, legendfontsize=12) 
    all_p = plot(p..., layout = (3, 2), size = (800, 800), linewidth = 2, thickness_scaling = 1.1, fontfamily = "Computer Modern")
    return all_p
end

function irfs(sr, lr, m_par; irf_horizon = 50)
    x = zeros(size(lr.LOMstate, 1), 1)
    x[sr.indexes_r.Z] = 100 * m_par.σ_Z

    MX = [I; lr.State2Control]
    irf_horizon = 40
    
    IRF = zeros(sr.n_par.ntotal_r, irf_horizon)

    for t = 1:irf_horizon
        IRF[:, t] = (MX * x)
        x = lr.LOMstate * x
    end

    return IRF
end