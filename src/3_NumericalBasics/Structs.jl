@doc raw"""
@make_struct_aggr(struct_name) 

Make `struct` `struct_name` with two fields for every variable name in `aggr_names`
(for steady state value and for deviation from it).

# Requires
(module) global `aggr_names`
"""
macro make_struct_aggr(struct_name)
	a_names 			= Symbol.(aggr_names)
	n_states 			= length(a_names)

	fields_states 	  	= [:($(a_names[i])::Int) for i = 1:n_states]
	fieldsSS_states   	= [:($(Symbol(a_names[i], "SS")) ::Int) for i = 1:n_states]

	esc(quote 
			struct $struct_name
					$(fieldsSS_states...)
					$(fields_states...)
			end
	end
	)
end

@doc raw"""
@make_struct(struct_name)

Make `struct` `struct_name` with two fields for every variable name in `s_names` (state variables)
and `c_names` (control variables), together with fields for distribution-states
and marginal value function-controls.

# Requires
(module) globals `state_names`, `control_names`
"""
macro make_struct(struct_name)
	# fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
	s_names 		    = Symbol.(state_names)
	n_states 			= length(s_names)
	c_names 			= Symbol.(control_names)
	n_controls 			= length(c_names)

	fields_states 	 	= [:($(s_names[i])::Int) for i = 1:n_states]
	fieldsSS_states   	= [:($(Symbol(s_names[i], "SS")) ::Int) for i = 1:n_states]
	fields_controls   	= [:($(c_names[i])::Int) for i = 1:n_controls]
	fieldsSS_controls 	= [:($(Symbol(c_names[i], "SS")) ::Int) for i =1:n_controls]

	# fields=[:($(Symbol(i, "::Int"))) for i in var_names]
	# fieldsSS=[:($(Symbol(i, "SS::Int"))) for i in var_names]
	esc(quote 
			struct $struct_name
				distr_m_SS	::Array{Int,1}
				distr_k_SS	::Array{Int,1}
				distr_y_SS	::Array{Int,1}
				COPSS     	::Array{Int,1}
				$(fieldsSS_states...)
				VmSS     	::Array{Int,1}
				VkSS     	::Array{Int,1}
				$(fieldsSS_controls...)
				distr_m		::Array{Int,1}
				distr_k		::Array{Int,1}
				distr_y		::Array{Int,1}
				COP    		::Array{Int,1}
				$(fields_states...)
				Vm     		::Array{Int,1}
				Vk     		::Array{Int,1}
				$(fields_controls...)
			end
		end
	)
end
struct SteadyResults
	XSS
	XSSaggr
	indexes
	indexes_aggr
	compressionIndexes
	n_par
	m_par
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
	nk
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
  