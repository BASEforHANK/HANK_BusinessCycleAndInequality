#------------------------------------------------------------------------------
# Header: load module
#------------------------------------------------------------------------------
# make sure that your pwd is set to the folder containing script and HANKEstim
# otherwise adjust the load path
# cd("HANK_BusinessCycleAndInequality/src")

# pre-process user inputs for model setup
include("3_NumericalBasics/PreprocessInputs.jl")
using BenchmarkTools, LinearAlgebra

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
e_set = HANKEstim.e_set;
# alternatively, load estimated parameters by running, e.g.,
# @load HANKEstim.e_set.save_posterior_file par_final e_set
# !! This file is not provided, because it is too large !!
# m_par = HANKEstim.Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(e_set.meas_error_input)])

################################################################################
# Comment in the following block to be able to go straight to plotting
################################################################################
# @load "7_Saves/steadystate.jld2" sr_full
# @load "7_Saves/linearresults.jld2" lr_full
# @load "7_Saves/reduction.jld2" sr_reduc lr_reduc
# @load HANKEstim.e_set.save_mode_file sr_mode lr_mode er_mode m_par_mode smoother_output
# @set! e_set.estimate_model = false

# Calculate Steady State at prior mode to find further compressed representation of Vm, Vk
sr_full = compute_steadystate(m_par)
jldsave("7_Saves/steadystate.jld2", true; sr_full) # true enables compression
# @load "7_Saves/steadystate.jld2" sr_full

lr_full = linearize_full_model(sr_full, m_par)
jldsave("7_Saves/linearresults.jld2", true; lr_full)
# @load "7_Saves/linearresults.jld2" lr_full

# Find sparse state-space representation
sr_reduc = model_reduction(sr_full, lr_full, m_par);
lr_reduc = update_model(sr_reduc, lr_full, m_par)
jldsave("7_Saves/reduction.jld2", true; sr_reduc, lr_reduc)
# @load "7_Saves/reduction.jld2" sr_reduc lr_reduc

println("One model solution takes")
#@load "7_Saves/linearresults.jld2" lr_full
#@load "7_Saves/reduction.jld2" sr_reduc lr_reduc
@set! sr_reduc.n_par.verbose = false
@btime lr_reduc = update_model(sr_reduc, lr_full, m_par)
@set! sr_reduc.n_par.verbose = true;

if e_set.estimate_model == true

        # warning: estimation might take a long time!
        er_mode, posterior_mode, smoother_mode, sr_mode, lr_mode, m_par_mode =
                find_mode(sr_reduc, lr_reduc, m_par)

        # Stores results in file e_set.save_mode_file 
        jldsave(HANKEstim.e_set.save_mode_file, true;
                posterior_mode, smoother_mode, sr_mode, lr_mode, er_mode, m_par_mode, e_set)
        # !! warning: the provided mode file does not contain smoothed covars (smoother_mode[4] and [5])!!
        # @load HANKEstim.e_set.save_mode_file posterior_mode sr_mode lr_mode er_mode m_par_mode smoother_mode e_set

        sr_mc, lr_mc, er_mc, m_par_mc, draws_raw, posterior, accept_rate,
        par_final, hessian_sym, smoother_output = montecarlo(sr_mode, lr_mode, er_mode, m_par_mode)

        # Stores results in file e_set.save_posterior_file 
        jldsave(HANKEstim.e_set.save_posterior_file, true;
                sr_mc, lr_mc, er_mc, m_par_mc, draws_raw, posterior, accept_rate,
                par_final, hessian_sym, smoother_output, e_set)

        # !! The following file is not provided !!
        #  @load HANKEstim.e_set.save_posterior_file sr_mc lr_mc er_mc  m_par_mc draws_raw posterior accept_rate par_final hessian_sym smoother_output e_set
end

###############################################################################################
# Graphical Model Output, functions not integrated in package
###############################################################################################
using Plots, VegaLite, DataFrames, FileIO, StatsPlots, CategoricalArrays, Flatten, Statistics
include("8_PostEstimation/IRFs_and_Decomp.jl")
# calculate uncompressed solution for comparison, needs either computed or loaded mode
@set! sr_full.n_par.further_compress = false
sr_mode_f = model_reduction(sr_full, lr_mode, m_par_mode);
lr_mode_f = update_model(sr_mode_f, lr_mode, m_par_mode)

select_variables = [:Y, :C, :I, :G, :TOP10Wshare] # Variables to be plotted
model_names = ["mode full", "mode reduction"] # Displayed names of models to be compared
# enter here the models, as tupel of tupels (sr, lr, e_set, m_par), to be compared
models_tupel = ((sr_mode_f, lr_mode_f, HANKEstim.e_set, m_par_mode),
        (sr_mode, lr_mode, HANKEstim.e_set, m_par_mode))
timeline = collect(1954.75:0.25:2019.75)
select_vd_horizons = [4 16 100] # horizons for variance decompositions
#recessions_vec   = [1957.5,1958.25,1960.25,1961.0,1969.75,1970.75,1973.75,1975.0,1980.0,1980.5,1981.5,1982.75,1990.5,1991.0, 2001.0, 2001.75, 2007.75, 2009.25]

IRFs, VDs, SHOCKs = compute_irfs_vardecomp(models_tupel, select_variables) # Compute IRFs for all models in tupel, all variables in select_variables

IRFs_plot = plot_irfs(IRFs, SHOCKs, select_variables, 20, model_names, 2; savepdf = false) # display IRFs and export as pdf

DF_V_Decomp = plot_vds(VDs, select_vd_horizons, model_names, SHOCKs, select_variables; savepdf = false) # Export Variance Decompositions as DataFrames and Plot using VegaLite

Historical_contrib, DF_H_Decomp, HD_plot = compute_hist_decomp(sr_mode, lr_mode, e_set, m_par_mode,
        smoother_output, select_variables, timeline; savepdf = false) # Produce historical contributions as Array and Data Frame and plot p

# Compute credible intervals for variance decomposition based on posterior draws
# !! The following file is not provided !!
# @load HANKEstim.e_set.save_posterior_file sr_mc lr_mc er_mc  m_par_mc draws_raw posterior accept_rate par_final hessian_sym smoother_output e_set
# models = ((sr_mc, lr_mc, HANKEstim.e_set, m_par_mc, draws),) # needs to be adjusted if more than one model
# VDcredint = compute_vardecomp_bounds(models, select_variables, select_vd_horizons, model_names)