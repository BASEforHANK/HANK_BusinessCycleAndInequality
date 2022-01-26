#------------------------------------------------------------------------------
# Header: load module
#------------------------------------------------------------------------------
# make sure that your pwd is set to the folder containing script and HANKEstim
# otherwise adjust the load path
# cd("HANK_BusinessCycleAndInequality/src")

# pre-process user inputs for model setup
include("3_NumericalBasics/PreprocessInputs.jl")
using BenchmarkTools, LinearAlgebra, Plots

push!(LOAD_PATH, pwd())
using HANKEstim

# set BLAS threads to the number of Julia threads.
# prevents BLAS from grabbing all threads on a machine
BLAS.set_num_threads(Threads.nthreads())

#------------------------------------------------------------------------------
# initialize parameters to priors to select coefficients of DCTs of Vm, Vk
# that are retained 
#------------------------------------------------------------------------------
m_par = ModelParameters()
priors = collect(metaflatten(m_par, prior)) # model parameters
par_prior = mode.(priors)
m_par = HANKEstim.Flatten.reconstruct(m_par, par_prior)

# alternatively, load estimated parameters from HANKX+ estimation and update model parameters
# @load "7_Saves/parameter_example.jld2" par_final
# m_par = HANKEstim.Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(HANKEstim.e_set.meas_error_input)])

# Calculate Steady State at prior mode to find further compressed representation of Vm, Vk
 sr_full = compute_steadystate(m_par)
 @save "7_Saves/steadystate.jld2" sr_full
#@load "7_Saves/steadystate.jld2" sr_full

 lr_full = linearize_full_model(sr_full, m_par)
 @save "7_Saves/linearresults.jld2" lr_full
#@load "7_Saves/linearresults.jld2" lr_full

# Find sparse state-space representation
 sr_reduc    = model_reduction(sr_full, lr_full, m_par);
 lr_reduc    = update_model(sr_reduc, lr_full, m_par)
 @save "7_Saves/reduction.jld2" sr_reduc lr_reduc
# println("One model solution takes")
# @set! sr_reduc.n_par.verbose = false
# @btime lr_reduc    = update_model(sr_reduc, lr_full, m_par)
# @set! sr_reduc.n_par.verbose = true

@load "7_Saves/reduction.jld2" sr_reduc lr_reduc

reduction_diagnostics = true
if reduction_diagnostics
    # plot some irfs to tfp (z) shock for full and reduced model


    function irfs(sr, lr; irf_horizon = 50)
        x0 = zeros(size(lr.LOMstate, 1), 1)
        x0[sr.indexes_r.Z] = 100 * m_par.σ_Z

        MX = [I; lr.State2Control]
        irf_horizon = 40
        x = x0 * ones(1, irf_horizon + 1)
        IRF_state_sparse = zeros(sr.n_par.ntotal_r, irf_horizon)

        for t = 1:irf_horizon
            IRF_state_sparse[:, t] = (MX * x[:, t])'
            x[:, t+1] = lr.LOMstate * x[:, t]
        end

        return IRF_state_sparse
    end

    IRF_state_sparse_r = irfs(sr_reduc, lr_reduc)
    IRF_state_sparse_o = irfs(sr_full, lr_full)

    p = Vector{Plots.Plot{Plots.GRBackend}}(undef, 6)
    p[1] = plot(IRF_state_sparse_r[sr_reduc.indexes_r.Z, :], title = "TFP (percent)", label = "reduced model")
    p[1] = plot!(IRF_state_sparse_o[sr_full.indexes_r.Z, :], label = "full model")
    p[2] = plot(IRF_state_sparse_r[sr_reduc.indexes_r.I, :], title = "Investment (percent)", label = "reduced model")
    p[2] = plot!(IRF_state_sparse_o[sr_full.indexes_r.I, :], label = "full model")
    p[3] = plot(IRF_state_sparse_r[sr_reduc.indexes_r.Y, :], title = "GDP (percent)", label = "reduced model")
    p[3] = plot!(IRF_state_sparse_o[sr_full.indexes_r.Y, :], label = "full model")
    p[4] = plot(IRF_state_sparse_r[sr_reduc.indexes_r.C, :], title = "Consumption (percent)", label = "reduced model")
    p[4] = plot!(IRF_state_sparse_o[sr_full.indexes_r.C, :], label = "full model")
    p[5] = plot(IRF_state_sparse_r[sr_reduc.indexes_r.TOP10Wshare, :], title = "Top 10 wealth share (percent)", label = "reduced model")
    p[5] = plot!(IRF_state_sparse_o[sr_full.indexes_r.TOP10Wshare, :], label = "full model")
    p[6] = plot(0)
    all_p = plot(p..., layout = (3, 2), size = (800, 800), linewidth = 2, thickness_scaling = 1.1, fontfamily = "Computer Modern")
    display(all_p)
end

if HANKEstim.e_set.estimate_model == true
    # warning: estimation might take a long time!
    er_mode, posterior_mode, smoother_output, sr_mode, lr_mode, m_par, e_set = find_mode(sr_reduc, lr_reduc, m_par)
    @save e_set.save_mode_file posterior_mode smoother_output sr_mode lr_mode er_mode m_par e_set
    # @load HANKEstim.e_set.save_mode_file sr_mode lr_mode er_mode m_par e_set
    # loading the mode only works with a full mode save file not our provided file
    montecarlo(sr_mode, lr_mode, er_mode, m_par) # Stores results in file e_set.save_posterior_file
end