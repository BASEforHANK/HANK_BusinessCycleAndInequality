###############################################################################################
# Compute IRFs and variance decomposition for set of models and variables passed to function
###############################################################################################
function compute_irfs_vardecomp(models, select_variables)

    max_horizon = 1000
    n_models = length(models)
    IRFs = Array{Array{Float64}}(undef, n_models)
    VDs  = Array{Array{Float64}}(undef, n_models)
    SHOCKs = Array{Symbol}(undef, 0)
    for j = 1: n_models
        e_set = models[j][3]
        union!(SHOCKs, e_set.shock_names)
    end
    n_shocks = length(SHOCKs)
    for j = 1:n_models
        sr    = models[j][1]
        lr    = models[j][2]
        e_set = models[j][3]
        m_par = models[j][4]
        selector = []
        isstate = zeros(Bool, length(select_variables))
        iter = 1
        for i in select_variables
            if i in Symbol.(sr.state_names)
                isstate[iter] .= true
            end
            iter += 1
            try
                append!(selector, getfield(sr.indexes_r, i))
            catch
                append!(selector, sr.ntotal_r + 1)
            end
        end

        IRFs_aux = zeros(length(selector), max_horizon, n_shocks)
        shock_number = 0
        for i in SHOCKs
            x = zeros(size(lr.LOMstate, 1))
            shock_number += 1
            try
                x[getfield(sr.indexes_r, i)] = getfield(m_par, Symbol("σ_", i))
            catch
                println("model ", j, " has no shock ", string(i))
            end
            MX = [I; lr.State2Control; zeros(1, sr.n_par.nstates_r)]
            for t = 1:max_horizon
                IRFs_aux[:, t, shock_number] = MX[selector, :] * x
                x[:] = lr.LOMstate * x
            end
        end
        IRFs_aux[isstate, 1:end-1, :] .= IRFs_aux[isstate, 2:end, :] # IRFs for state variables represent end-of-period values
        # Dimensions of IRF: variable x time x shock
        # VARdecomp = zeros(size(IRFs_aux))
        
        VARdecomp = cumsum(IRFs_aux .^ 2.0, dims = 2) ./ (sum(cumsum(IRFs_aux .^ 2.0, dims = 2), dims = 3) .+ 10.0 * eps())
        
        VDs[j] = 100.0 .* VARdecomp
        IRFs[j] = 100.0 .* IRFs_aux

    end
    return IRFs, VDs, SHOCKs

end

function compute_vardecomp_bounds(models, select_variables, select_vd_horizons, model_names; n_replic = 1000, percentile_bounds = (0.05, 0.95))

    VDorig, SHOCKs = compute_irfs_vardecomp(models, select_variables)[2:3]
    max_horizon = maximum(select_vd_horizons)
    n_models = length(models)
    VD_lower = Array{Array{Float64}}(undef, n_models)
    VD_upper = Array{Array{Float64}}(undef, n_models)

    n_shocks = length(SHOCKs)
    for j = 1:n_models
        sr = models[j][1]
        lr = models[j][2]
        e_set = models[j][3]
        m_par = models[j][4]
        draws = models[j][5]
        selector = []
        isstate = zeros(Bool, length(select_variables))
        iter = 1
        for i in select_variables
            if i in Symbol.(sr.state_names)
                isstate[iter] .= true
            end
            iter += 1
            try
                append!(selector, getfield(sr.indexes_r, i))
            catch
                append!(selector, sr.ntotal_r + 1)
            end
        end
        VDaux = hcat([zeros(n_replic) for j = 1:length(select_variables), i = 1:max_horizon, k = 1:n_shocks])
        draw_ind = rand(1:size(draws, 1), n_replic)
        for s = 1:length(draw_ind)
            dd = draws[draw_ind[s], :]
            if e_set.me_treatment != :fixed
                m_par = Flatten.reconstruct(m_par, dd[1:length(par_final)-length(e_set.meas_error_input)])
            else
                m_par = Flatten.reconstruct(m_par, dd)
            end
            lr = update_model(sr, lr, m_par)
    
            IRFs_aux = zeros(length(selector), max_horizon, n_shocks)
            shock_number = 0
            for i in SHOCKs
                x = zeros(size(lr.LOMstate, 1))
                shock_number += 1
                try
                    x[getfield(sr.indexes_r, i)] = getfield(m_par, Symbol("σ_", i))
                catch
                    println("model ", j, " has no shock ", string(i))
                end
                MX = [I; lr.State2Control; zeros(1, sr.n_par.nstates_r)]
                for t = 1:max_horizon
                    IRFs_aux[:, t, shock_number] = MX[selector, :] * x
                    x[:] = lr.LOMstate * x
                end
            end
            IRFs_aux[isstate, 1:end-1, :] .= IRFs_aux[isstate, 2:end, :] # IRFs for state variables represent end-of-period values
            # Dimensions of IRF: variable x time x shock
    
            VDauxaux = 100 .* cumsum(IRFs_aux .^ 2.0, dims = 2) ./ (sum(cumsum(IRFs_aux .^ 2.0, dims = 2), dims = 3) .+ 10.0 * eps())
            for d1 = 1:length(select_variables)
                for d2 = 1:max_horizon
                    for d3 = 1:n_shocks
                        VDaux[d1, d2, d3][s] = VDauxaux[d1, d2, d3]
                    end
                end
            end
        end
    
        VD_lower[j] = quantile.(VDaux, percentile_bounds[1])
        VD_upper[j] = quantile.(VDaux, percentile_bounds[2])
    end
    n_total_entries = n_models * length(select_variables) * length(select_vd_horizons) * n_shocks
    modelname_vec = Vector{String}(undef, n_total_entries)
    variablename_vec = Vector{Symbol}(undef, n_total_entries)
    horizonname_vec = Vector{Int}(undef, n_total_entries)
    shockname_vec = Vector{Symbol}(undef, n_total_entries)
    VD_vec = Vector{Float64}(undef, n_total_entries)
    VDlow_vec = Vector{Float64}(undef, n_total_entries)
    VDup_vec = Vector{Float64}(undef, n_total_entries)
    count = 1
    for d0 = 1:n_models
        for d1 = 1:length(select_variables)
            for d2 in select_vd_horizons
                for d3 = 1:n_shocks
                    modelname_vec[count] = model_names[d0]
                    variablename_vec[count] = select_variables[d1]
                    horizonname_vec[count] = d2
                    shockname_vec[count] = SHOCKs[d3]
                    VD_vec[count] = copy(VDorig[d0][d1, d2, d3])
                    VDlow_vec[count] = copy(VD_lower[d0][d1, d2, d3])
                    VDup_vec[count] = copy(VD_upper[d0][d1, d2, d3])

                    count += 1
                end
            end
        end
    end
    VDdf = DataFrame(model = modelname_vec, variable = variablename_vec, horizon = horizonname_vec,
        shock = shockname_vec, lower_bound = VDlow_vec, point_estimate = VD_vec, upper_bound = VDup_vec)

    return VDdf

end

shock_color = ["#ff0000", "#ffcccc", "#ff9900", "#ffebcc", "#009900", "#ccffcc", "#F7DBFF", "#0000ff", "#ccccff"]
ds_color = ["#0000ff", "#ff0000"]

###############################################################################################
# Plot IRFs
###############################################################################################
function plot_irfs(IRFs, SHOCKs, select_variables, horizon, model_names::Array{String}, n_plotcol; savepdf = false)
    styles = [:solid :dash :dashdot :dashdotdot :dot :dash :solid]
    colorlist = [:black, :red, :blue, :green, :orange, :purple, :yellow]
    pvec = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(SHOCKs))
    counts = 0
    for s in SHOCKs
        counts += 1
        countv = 0
        p = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(select_variables) + 1)
        for v in select_variables
            countv += 1
            p[countv] = plot()
            countm = 0
            for irf in IRFs
                countm += 1
                p[countv] = plot!(irf[countv, 1:horizon, counts], linewidth = 2,
                    linestyle = styles[countm], palette = colorlist, legend = false, title = v)
            end
        end
        # Add a plot that only contains the legend
        p[end] = plot(Matrix{Missing}(undef, 3, length(model_names)), label = [string(z) for k = 1:1, z in model_names], legend = :inside,
            linestyle = styles, palette = colorlist,
            framestyle = :none, legendfontsize = 14, title = string("Response to: ", s))
        # Combine in a plot with sublots
        pvec[counts] = plot(p..., layout = (ceil(Int, length(p) / n_plotcol), n_plotcol), size = (800, 800),
            linewidth = 2, thickness_scaling = 1.1, fontfamily = "Computer Modern")
    end
    if savepdf
        for j = 1:length(SHOCKs)
            pvec[j] |> save(string("8_PostEstimation/Figures/IRFs_to_", SHOCKs[j], ".pdf"))
        end
    end
    display.(pvec)
    return pvec
end

###############################################################################################
# Plot variance decomposition
###############################################################################################
function plot_vds(VDs, select_vd_horizons, model_names, SHOCKs, select_variables; savepdf = false)

    SHOCKs = CategoricalArray(string.(SHOCKs))
    levels!(SHOCKs,SHOCKs)
    select_variables = CategoricalArray(string.(select_variables))
    levels!(select_variables,select_variables)
    model_names = CategoricalArray(model_names)
    levels!(model_names,model_names)
    # construct DataFrame that stacks variance decomposition
    df = vcat([HANKEstim.DataFrame(Model = model_names[i], Horizon = h,
        VarianceDecomposition = VDs[i][j, h, k], Shock = SHOCKs[k], Variable = select_variables[j]) 
        for i = 1:length(model_names), h in select_vd_horizons, k = 1:length(SHOCKs), j = 1:length(select_variables)]...)

    # plot variance decomposition
    for j in select_variables
        p = ((df[df.Variable .== j, :] |> @vlplot(:bar,
            x = {:VarianceDecomposition,
                title = nothing,
                scale = {domain = [0, 100], nice = false}},
            color = {:Shock, scale = {range = shock_color}, sort = SHOCKs},
            row = {"Model:n", title = "Model", sort = model_names},
            column = {"Horizon:n", title = string("Variance Decomposition for ", j, ", Forecast Horizon:"), sort = select_vd_horizons},
            background = :white,
            order = "siteOrder:o",
            width = 200,
            height = 50)))
        if savepdf
            p |> save(string("8_PostEstimation/Figures/CVD_of_", j, ".pdf"))
        end
        display(p)
    end

    return df
end

###############################################################################################
# Plot historical decomposition
###############################################################################################
function compute_hist_decomp(sr, lr, e_set, m_par, smoother_output, select_variables,timeline; savepdf=false)
    
    SHOCKs = e_set.shock_names
    n_shocks = length(SHOCKs)
    n_select_var = length(select_variables)
    T = size(smoother_output[6], 2)
    IRFs = Array{Float64}(undef, n_select_var, T, n_shocks + 1)
    
    # select variables of interest from all variables
    selector = []
    isstate = zeros(Bool, length(select_variables))
    iter = 1
    for i in select_variables
        if i in Symbol.(sr.state_names)
            isstate[iter] .= true
        end
        iter += 1
        try
            append!(selector, getfield(sr.indexes_r, i))
        catch
            append!(selector, sr.n_par.ntotal_r + 1)
        end
    end

    MX = [I; lr.State2Control]

    # effect of initial condition
    IRFs_aux = zeros(length(selector), T)
    x = smoother_output[3][:,1] 
    for t = 1:T
        IRFs_aux[:, t] = MX[selector, :] * x
        x[:] = lr.LOMstate * x
    end
    IRFs[:,:,n_shocks+1] = IRFs_aux

    # loop through shocks, calculate effect of each shock
    for j  = 1:n_shocks
        IRFs_aux = zeros(length(selector), T)
        x = zeros(size(lr.LOMstate, 1))
        i = SHOCKs[j] 
        shock_index = getfield(sr.indexes_r, i)
        for t = 1:T
            IRFs_aux[:, t] = MX[selector, :] * x
            x[:] = lr.LOMstate * x
            x[shock_index] += smoother_output[6][shock_index,t] # shock in "t" moves observables in "t+1" 
        end
        IRFs[:,:,j] = IRFs_aux
    end

    SHOCKsPlus = push!(copy(SHOCKs), :initial)
    SHOCKsPlus = CategoricalArray(string.(SHOCKsPlus))
    levels!(SHOCKsPlus, SHOCKsPlus)
    select_variables = CategoricalArray(string.(select_variables))
    levels!(select_variables ,select_variables )
    HistDecDF  = vcat([DataFrame(Time = timeline[t], 
                                Variable = select_variables[v], 
                                Contribution = IRFs[v,t,s], 
                                Shock = SHOCKsPlus[s]) for t=1:T, v=1:n_select_var, s=1:n_shocks+1]...)
    
    colorlist = ["#ff0000", "#ffcccc", "#ff9900", "#ffebcc", "#009900", "#ccffcc", "#F7DBFF", "#0000ff", "#ccccff", :grey]

    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(select_variables) + 1)
    # rec = vspan(recessions_vec, fill=:grey85, linecolor=:grey85, label="")
    countj=0
    for j in select_variables
        countj+=1

        p[countj] = @df HistDecDF[HistDecDF.Variable .== j,:] groupedbar(:Time, :Contribution, group = :Shock, linecolor =false, title = j,bar_position = :stack, legend=false, palette = colorlist)
        p[countj] = @df combine(groupby(HistDecDF[HistDecDF.Variable .== j,:], :Time), :Contribution => sum) plot!(:Time, :Contribution_sum, label = string(j), linewidth =1, linecolor=:black, legend=false)
    end
    p[end] = bar(Matrix{Missing}(undef, 3, n_shocks+1), label = [string(z) for k = 1:1, z in SHOCKsPlus],linecolor =false, legend = :inside,
                    palette = colorlist, framestyle = :none, legendfontsize = 10)
    p[end] = plot!(Matrix{Missing}(undef, 3, 1), label = "Total", legend = :inside,
                    linecolor=:black, framestyle = :none, legendfontsize = 10)
    if savepdf
        plot(p...) |> save(string("8_PostEstimation/Figures/HistD.pdf"))
    end
    display(plot(p...))
    return IRFs, HistDecDF, p
end