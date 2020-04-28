# make sure that your pwd is set to the folder containing script and HANKEstim
# otherwise adjust the load path

push!(LOAD_PATH, pwd())
using HANKEstim, JLD2

# sr = compute_steadystate()
# save_steadystate(sr)
sr = load_steadystate()

# lr = linearize_full_model(sr)
# @save "linearresults.jld2" lr
@load "Saves/linearresults.jld2"

# warning: estimation might take a long time!
# er = find_mode(sr,lr)
# @save "estimresults.jld2" er
# @load "estimresults.jld2"
# montecarlo(sr,lr,er)

# plot some irfs to tfp (z) shock
using LinearAlgebra, Plots
x0      = zeros(size(lr.LOMstate,1), 1)
x0[sr.indexes.Z] = sr.m_par.σ_Z

MX = [I; lr.State2Control]
nlag=20
x = x0*ones(1,nlag+1)
IRF_state_sparse = zeros(sr.indexes.profits, nlag)
for t = 1:nlag
        IRF_state_sparse[:, t]= (MX*x[:,t])'
        x[:,t+1] = lr.LOMstate * x[:,t]
end
plt1 = plot(IRF_state_sparse[sr.indexes.Z,:],label="TFP (percent)", reuse=false)
plt1 =plot!(IRF_state_sparse[sr.indexes.I,:],label="Investment (percent)")
plt1 =plot!(IRF_state_sparse[sr.indexes.Y,:],label="Output (percent)")
plt1 =plot!(IRF_state_sparse[sr.indexes.C,:],label="Consumption (percent)")
