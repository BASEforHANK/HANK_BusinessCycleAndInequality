# make sure that your pwd is set to the folder containing script and HANKEstim
# otherwise adjust the load path

push!(LOAD_PATH, pwd())
using HANKEstim

#initialize model parameters
m_par = ModelParameters()
# load estimated parameters from HANKX+ estimation and update model parameters
HANKEstim.@load "Saves/HANKXplus_postmean.jld2" par_final
m_par = HANKEstim.Flatten.reconstruct(m_par, par_final[1:length(par_final) - length(HANKEstim.e_set.meas_error_input)])

# Calculate Steady State
sr    = compute_steadystate(m_par)
HANKEstim.@save "Saves/steadystate.jld2" sr
# HANKEstim.@load "Saves/steadystate.jld2" sr

lr = linearize_full_model(sr, m_par)
HANKEstim.@save "Saves/linearresults.jld2" lr
# HANKEstim.@load "Saves/linearresults.jld2"

# plot some irfs to tfp (z) shock
using LinearAlgebra, Plots
x0      = zeros(size(lr.LOMstate,1), 1)
x0[sr.indexes.Z] = 100 * m_par.σ_Z
# x0[sr.indexes.σ] = 100 * m_par.σ_Sshock

MX = [I; lr.State2Control]
irf_horizon = 40
x = x0 * ones(1, irf_horizon + 1)
IRF_state_sparse = zeros(sr.n_par.ntotal, irf_horizon)
for t = 1:irf_horizon
        IRF_state_sparse[:, t]= (MX * x[:, t])'
        x[:, t+1] = lr.LOMstate * x[:, t]
end
plt1 = plot(IRF_state_sparse[sr.indexes.Z,:],label="TFP (percent)", reuse=false)
plt1 = plot!(IRF_state_sparse[sr.indexes.I,:],label="Investment (percent)")
plt1 = plot!(IRF_state_sparse[sr.indexes.Y,:],label="Output (percent)")
plt1 = plot!(IRF_state_sparse[sr.indexes.C,:],label="Consumption (percent)")

if HANKEstim.e_set.estimate_model == true
    # warning: estimation might take a long time!
    er = find_mode(sr, lr, m_par)
    # loading the mode only works with a full mode save file not our provided file
    # er = load_mode(sr; file = HANKEstim.e_set.save_mode_file)

    montecarlo(sr, lr, er, m_par)
end