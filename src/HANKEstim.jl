__precompile__(false)
# Code runs on Julia 1.3
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------
# Packages used: Plots Distributions BenchmarkTools JLD2 FileIO DataFrames ForwardDiff
# SparseArrays LinearAlgebra Random LaTeXStrings Arpack SpecialFunctions FFTW
# Parameters Setfield MCMCChains StatsPlots Optim CSV OrderedCollections
# Flatten FieldMetadata JSON

module HANKEstim

using Plots, Distributions, BenchmarkTools, JLD2, FileIO, DataFrames, ForwardDiff
using SparseArrays, LinearAlgebra, Random, LaTeXStrings
using Arpack: eigs; using SpecialFunctions: erf; using FFTW: dct
using Parameters, Setfield, MCMCChains, StatsPlots, Optim, CSV, OrderedCollections
using Flatten; import Flatten: flattenable
using FieldMetadata; import FieldMetadata: prior, label
using JSON

export LinearResults, linearize_full_model, EstimResults, find_mode, load_mode, montecarlo,
        find_SS, @writeXSS, @make_fn, @make_fnaggr, @make_struct, @make_struct_aggr,
        mylinearinterpolate3, Ksupply, Fsys, SGU, Fsys_agg, SGU_estim, mode_finding,
        likeli, kalman_filter, kalman_filter_smoother, measurement_error, rwmh,
        ModelParameters, NumericalParameters, EstimationSettings,
        load_steadystate, save_steadystate, compute_steadystate, SteadyResults,
        Tauchen, EGM_policyupdate, Kdiff, distrSummaries, @generate_equations,
        @make_deriv, @make_deriv_estim, prioreval

include("input_aggregate_names.jl")

aggr_names = [state_names; control_names]

distr_names=["GiniW", "GiniC", "GiniX", "GiniI", "sdlgC", "P9010C", "I90share",
"I90sharenet", "P9010I", "w90share", "P10C", "P50C", "P90C"]

n_FD = 5 # number of derivatives to be calculated simultaneously,
# optimal choice depends on CPU and memory of machine


# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("Structs.jl")
include("Estimation/prior.jl")

e_set = EstimationSettings()
@make_struct IndexStruct state_names control_names
@make_struct_aggr IndexStructAggr aggr_names

include("NumericalBasics.jl")
include("HetAgentsFcns.jl")
include("LinearizationFunctions.jl")
include("Estimation.jl")

@make_fn produce_indexes state_names control_names
@make_fnaggr produce_indexes_aggr aggr_names

struct SteadyResults
  XSS
  XSSaggr
  indexes
  indexes_aggr
  compressionIndexes
  Copula
  n_par
  m_par
  d
  CDF_SS
  CDF_m
  CDF_k
  CDF_y
  distrSS
end

struct LinearResults
  State2Control
  LOMstate
  A
  B
  SolutionError
end

struct EstimResults
  par_final
  hessian_final
  meas_error
  meas_error_std
  parnames
  Data
  Data_missing
  H_sel
  priors
end

@doc raw"""
    save_steadystate(XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, m_par, d,
    CDF_SS, CDF_m, CDF_k, CDF_y, distrSS)

Save global variables of steady state to files `steadystate.jld2`,`m_par.json` and `n_par.json`.
"""
function save_steadystate(sr::SteadyResults)
  save("Saves/steadystate.jld2",Dict(
                                "XSS" => sr.XSS,
                                "XSSaggr" => sr.XSSaggr,
                                "compressionIndexes" => sr.compressionIndexes,
                                "d" => sr.d,
                                "CDF_SS" => sr.CDF_SS,
                                "CDF_m" => sr.CDF_m,
                                "CDF_k" => sr.CDF_k,
                                "CDF_y" => sr.CDF_y,
                                "distrSS" => sr.distrSS
    ))
  structs = [:m_par,:n_par]
  for istr = 1:2
    filestr = string("Saves/", String(structs[istr]),".json")
    if isfile(filestr)
      rm(filestr)
    end
    open(filestr,"w") do f
      JSON.print(f,getfield(sr,structs[istr]),4)
    end
  end
end

@doc raw"""
    load_steadystate()

Load global variables of steady state from files (filenames as in [`save_steadystate()`](@ref)).

# Returns
same returns as [`find_SS()`](@ref)
"""
function load_steadystate()
  dictmpar = JSON.parsefile("Saves/m_par.json")
  mpairs = Array{Pair{Symbol,Any},1}()
  for (k,v) in dictmpar
    push!(mpairs,Pair{Symbol,Any}(Symbol(k),v))
  end
  m_par = ModelParameters(;mpairs...)
  dictnpar = JSON.parsefile("Saves/n_par.json")
  npairs = Array{Pair{Symbol,Any},1}()
  for (k,v) in dictnpar
    # Deal with Symbol-values (that are saved as Strings) and Array-values separately.
    if String == typeof(v) || Array{Any,1} == typeof(v)
      continue
    end
    push!(npairs,Pair{Symbol,Any}(Symbol(k),v))
  end
  n_par = NumericalParameters(;grid_y = dictnpar["grid_y"],
                               grid_k = dictnpar["grid_k"],
                               grid_m = dictnpar["grid_m"],
                               sol_algo = Symbol(dictnpar["sol_algo"]),
                               aggr_names = dictnpar["aggr_names"],
                               bounds_y = dictnpar["bounds_y"],
                               npairs...
                             )
  # Deal with all complicated values manually.
  # Infer size of multi-dim. arrays from other parameters.
  @set! n_par.Π = reshape(vcat(dictnpar["Π"]...),(n_par.ny,n_par.ny))
  @set! n_par.mesh_y = reshape(vcat(vcat(dictnpar["mesh_y"]...)...),
      (n_par.nm,n_par.nk,n_par.ny))
  @set! n_par.mesh_k = reshape(vcat(vcat(dictnpar["mesh_k"]...)...),
      (n_par.nm,n_par.nk,n_par.ny))
  @set! n_par.mesh_m = reshape(vcat(vcat(dictnpar["mesh_m"]...)...),
      (n_par.nm,n_par.nk,n_par.ny))
  @set! n_par.dist_guess = reshape(vcat(vcat(dictnpar["dist_guess"]...)...),
      (n_par.nm,n_par.nk,n_par.ny))

  @load "Saves/steadystate.jld2" XSS XSSaggr compressionIndexes d CDF_SS CDF_m CDF_k CDF_y distrSS
  # produce indexes to access XSS etc.
  indexes = produce_indexes(n_par, compressionIndexes[1], compressionIndexes[2])
  indexes_aggr = produce_indexes_aggr(n_par)
  Copula(x::Vector,y::Vector,z::Vector) = mylinearinterpolate3(CDF_m, CDF_k, CDF_y,
                                                              CDF_SS, x, y, z)
  return SteadyResults(XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, m_par, d,
  CDF_SS, CDF_m, CDF_k, CDF_y, distrSS)
end

@doc raw"""
    compute_steadystate()

Compute steady state.

# Returns
`struct` `SteadyResults`, containing returns of [`find_SS()`](@ref)
"""
function compute_steadystate()
  XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, m_par, d,
  CDF_SS, CDF_m, CDF_k, CDF_y, distrSS  = find_SS()
  # Convention: profits is the last control in the list of control variables
  ntotal                  = indexes.profits
  @set! n_par.ntotal      = ntotal
  @set! n_par.ncontrols   = length(compressionIndexes[1]) + length(compressionIndexes[2]) + n_par.naggrcontrols

  @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
  @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
  return SteadyResults(XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, Copula, n_par, m_par, d,
  CDF_SS, CDF_m, CDF_k, CDF_y, distrSS)
end

@doc raw"""
    linearize_full_model()

Linearize the full solution (i.e. with idiosyncratic states and controls) around the steady state,
using [`SGU()`](@ref).

# Returns
`struct` `LinearResults`, containing
- `A::Array{Float64,2}`,`B::Array{Float64,2}`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`]
    and `XPrime` [`A`]
- `State2Control::Array{Float64,2}`: observation equation
- `LOMstate::Array{Float64,2}`: state transition equation
"""
function linearize_full_model(sr::SteadyResults)
    A=zeros(sr.n_par.ntotal,sr.n_par.ntotal)
    B=zeros(sr.n_par.ntotal,sr.n_par.ntotal)
    State2Control, LOMstate, SolutionError, nk, A, B = SGU(sr.XSS, copy(A), copy(B), sr.m_par, sr.n_par, sr.indexes, sr.Copula, sr.compressionIndexes, sr.distrSS; estim = false)
    return LinearResults(State2Control, LOMstate, A, B, SolutionError)
end

@doc raw"""
    find_mode(sr,mr)

Find parameter that maximizes likelihood of data given linearized model `mr`.

# Arguments
- `sr::SteadyResults`
- `lr::LinearResults`

# Returns
`struct` `EstimResults`, containing all returns of [`mode_finding()`](@ref)
"""
function find_mode(sr::SteadyResults,lr::LinearResults)
    par_final, hessian_final, posterior_mode, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors =
        mode_finding(sr.XSSaggr, lr.A, lr.B, sr.indexes, sr.indexes_aggr, sr.Copula,sr.distrSS, sr.compressionIndexes, sr.m_par, sr.n_par, e_set)
    return EstimResults(par_final, hessian_final, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors)
end

function load_mode(sr::SteadyResults;file::String = e_set.mode_start_file)
    @load file par_final hessian_final hessian_sym meas_error meas_error_std parnames Data Data_missing H_sel priors

    # Load data
    Data_temp = CSV.read(e_set.data_file; missingstring="NaN")
    data_names_temp = names(Data_temp)
    for i in data_names_temp
      name_temp = get(e_set.data_rename, i, :none)
      if name_temp == :none

      else
        rename!(Data_temp, i=>name_temp)
      end
    end

    # reorder variables so that first-differences come first
    observed_vars = [e_set.observed_vars_input[e_set.growth_rate_select]; e_set.observed_vars_input[.!e_set.growth_rate_select]]
    Data = Matrix(Data_temp[observed_vars])
    Data_missing = ismissing.(Data)
    nobservables = size(Data)[2]

    n_growth_rates = sum(e_set.growth_rate_select) # number of growth rate variables
    if e_set.fd_flag == true
      H_sel = zeros(nobservables, 2 * (sr.n_par.nstates + sr.n_par.ncontrols))
      for i in eachindex(observed_vars)
        if i <= n_growth_rates
          H_sel[i,  getfield(indexes_,(observed_vars[i]))]  = 1.0
        else
          H_sel[i, sr.n_par.nstates + sr.n_par.ncontrols + getfield(sr.indexes,(observed_vars[i]))]  = 1.0
        end
      end
    else
      H_sel = zeros(nobservables, sr.n_par.nstates + sr.n_par.ncontrols)
      for i in eachindex(observed_vars)
        H_sel[i,   getfield(sr.indexes,(observed_vars[i]))]  = 1.0
      end
    end
    return EstimResults(par_final, hessian_final, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors)
end

@doc raw"""
    montecarlo(mr,er;file=e_set.save_posterior_file)

Sample posterior of parameter vector with [`rwmh()`](@ref), take sample mean as
parameter estimate, and save all results in `file`.

# Arguments
- `sr::SteadyResults`
- `mr::LinearResults`
- `er::EstimResults`
"""
function montecarlo(sr::SteadyResults,lr::LinearResults, er::EstimResults;file::String = e_set.save_posterior_file)
    hessian_sym = Symmetric(nearest_spd(inv(er.hessian_final)))

    # init_draw, init_success, init_iter = multi_chain_init(par_final, hessian_sym, n_par, Data, Data_missing,
    #  H_sel, XSSaggr,A,B, indexes, indexes_aggr,	m_par,  e_set, Copula, distrSS, compressionIndexes,
    #  priors, meas_error, meas_error_std)
    #  par_final = init_draw
    # if init_success == false
    #   error("Couldn't find initial value that produces posterior")
    # end
    # @set! e_set.save_posterior_file = "../Saves/hank_2asset_mcmc_1001_baseline_ineq_chain4.jld2"

    draws_raw, posterior, accept_rate = rwmh(er.par_final, hessian_sym, sr.n_par, er.Data, er.Data_missing,
    er.H_sel, sr.XSSaggr,lr.A, lr.B, lr.LOMstate, lr.State2Control, sr.indexes,sr.indexes_aggr, sr.m_par, e_set,
    sr.Copula, sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std)

    draws = draws_raw[e_set.burnin+1:end, :]

    ##
    parnames_ascii = collect(metaflatten(sr.m_par, label))
    if e_set.me_treatment != :fixed
    for i in eachindex(e_set.meas_error_input)
        push!(parnames_ascii, string("sigma_me_", e_set.meas_error_input[i]))
    end
    end

    # need to ] add MCMCChains#master to get the sorted keyword
    chn = Chains(reshape(draws, (size(draws)...,1)), [string(parnames_ascii[i]) for i = 1:length(parnames_ascii)], sorted=false)
    chn_summary = summarize(chn)
    par_final = chn_summary[:,:mean];

    ##
    if e_set.me_treatment != :fixed
    m_par = Flatten.reconstruct(sr.m_par, par_final[1:length(par_final) - length(er.meas_error)])
    else
    m_par = Flatten.reconstruct(sr.m_par, par_final)
    end
    State2Control, LOMstate, alarm_sgu, nk = SGU_estim(sr.XSSaggr,  lr.A, lr.B, m_par, sr.n_par, sr.indexes, sr.indexes_aggr, sr.distrSS; estim = true)

    smoother_output = likeli(par_final, er.Data, er.Data_missing, er.H_sel, sr.XSSaggr, lr.A, lr.B, sr.indexes, sr.indexes_aggr, m_par, sr.n_par, e_set,
    sr.Copula, sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std; smoother = true)
    ##
    save(file,Dict(
                    "draws_raw" => draws_raw,
                    "accept_rate" => accept_rate,
                    "aggr_names_ascii" => aggr_names_ascii,
                    "aggr_names" => aggr_names,
                    "state_names_ascii" => state_names_ascii,
                    "state_names" => state_names,
                    "control_names" => control_names,
                    "control_names_ascii" => control_names_ascii,
                    "e_set" => e_set,
                    "Data" => er.Data,
                    "State2Control" => State2Control,
                    "LOMState" => LOMstate,
                    "parnames" => er.parnames,
                    "indexes" => sr.indexes,
                    "indexes_aggr" => sr.indexes_aggr,
                    "par_final" => par_final,
                    "m_par" => m_par,
                    "A" => lr.A,
                    "B" => lr.B,
                    "hessian_sym" => hessian_sym,
                    "n_par" => sr.n_par,
                    "draws" => draws,
                    "posterior" => posterior,
                    "hessian_final" => er.hessian_final,
                    "meas_error" => er.meas_error,
                    "meas_error_std" => er.meas_error_std,
                    "Data_missing" => er.Data_missing,
                    "H_sel" => er.H_sel,
                    "priors" => er.priors,
                    "smoother_output" => smoother_output
    ))
end

end # module HANKEstim
