function incomes(n_par, m_par, mcw, A, q, RB, τprog, τlev, H, Ht, π,r,w,N, profits, unionprofits, av_tax_rate)


   # profits        = (1.0 -mc) .* Y
    tax_prog_scale = (m_par.γ + m_par.τ_prog)/((m_par.γ + τprog))
   # unionprofits   = w * N * (1.0 - mcw)
    eff_int      = ((RB .* A) ./ π        .+ (m_par.Rbar .* (n_par.mesh_m.<=0.0)))  # effective rate (need to check timing below and inflation)
    
    GHHFA = ((m_par.γ + τprog)/(m_par.γ+1)) # transformation (scaling) for composite good
    inc   = [
                GHHFA.*τlev.*((n_par.mesh_y/H).^tax_prog_scale .*mcw.*w.*N./(Ht)).^(1.0-τprog).+
                (unionprofits).*(1.0 .- av_tax_rate).* n_par.HW, # incomes of workers adjusted for disutility of labor
                (r .- 1.0).* n_par.mesh_k, # rental income
                eff_int .* n_par.mesh_m, # liquid asset Income
                n_par.mesh_k .* q,
                τlev.*(mcw.*w.*N.*n_par.mesh_y./ H).^(1.0-τprog).*((1.0 - τprog)/(m_par.γ+1)),
                τlev.*((n_par.mesh_y/H).^tax_prog_scale .*mcw.*w.*N./(Ht)).^(1.0-τprog) .+ 
                unionprofits.*(1.0 .- av_tax_rate).* n_par.HW# capital liquidation Income (q=1 in steady state)
                ] 
    inc[1][:,:,end].= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes
    inc[5][:,:,end].= 0.0
    inc[6][:,:,end].= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes

    incgross = [
                ((n_par.mesh_y/H).^tax_prog_scale .*mcw.*w.*N./(Ht)).+
                (unionprofits).* n_par.HW,
                (r .- 1.0).* n_par.mesh_k,                                      # rental income
                eff_int .* n_par.mesh_m,                                        # liquid asset Income
                n_par.mesh_k .* q,
                ((n_par.mesh_y/H).^tax_prog_scale .*mcw.*w.*N./(Ht))            # capital liquidation Income (q=1 in steady state)
                ]           
    incgross[1][:,:,end].= (n_par.mesh_y[:,:,end] .* profits)
    incgross[5][:,:,end].= (n_par.mesh_y[:,:,end] .* profits)

   
    # taxrev          = (incgross[5]-inc[6]) # tax revenues w/o tax on union profits
    # incgrossaux     = incgross[5][1,1,:]
    # av_tax_rate     = dot(distr_y, taxrev[1,1,:])./(dot(distr_y,incgrossaux))
    
    # inc[6]          = τlev.*((n_par.mesh_y/H).^tax_prog_scale .*mcw.*w.*N./(Ht)).^(1.0-τprog) .+ 
    #                   unionprofits.*(1.0 .- av_tax_rate).* n_par.HW
    # inc[6][:,:,end].= τlev.*(n_par.mesh_y[:,:,end] .* profits).^(1.0-τprog) # profit income net of taxes

return incgross, inc, eff_int
end