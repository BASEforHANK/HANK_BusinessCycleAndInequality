@doc raw"""
@make_struct_aggr(struct_name,a_names) 

Make `struct` `struct_name` with two fields for every variable name in `a_names`
(for steady state value and for deviation from it).
"""
macro make_struct_aggr(struct_name, a_names)
	a_names 			= Symbol.(eval((a_names)))
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
@make_struct(struct_name,s_names,c_names)

Make `struct` `struct_name` with two fields for every variable name in `s_names` (state variables)
and `c_names` (control variables), together with fields for distribution-states
and marginal value function-controls.
"""
macro make_struct(struct_name, s_names, c_names)
	# fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
	state_names 		= Symbol.(eval((s_names)))
	n_states 			= length(state_names)
	control_names 		= Symbol.(eval((c_names)))
	n_controls 			= length(control_names)

	fields_states 	 	= [:($(state_names[i])::Int) for i = 1:n_states]
	fieldsSS_states   	= [:($(Symbol(state_names[i], "SS")) ::Int) for i = 1:n_states]
	fields_controls   	= [:($(control_names[i])::Int) for i = 1:n_controls]
	fieldsSS_controls 	= [:($(Symbol(control_names[i], "SS")) ::Int) for i =1:n_controls]

	# fields=[:($(Symbol(i, "::Int"))) for i in var_names]
	# fieldsSS=[:($(Symbol(i, "SS::Int"))) for i in var_names]
	esc(quote 
			struct $struct_name
				distr_m_SS	::Array{Int,1}
				distr_k_SS	::Array{Int,1}
				distr_y_SS	::Array{Int,1}
				DSS     	::Array{Int,1}
				$(fieldsSS_states...)
				VmSS     	::Array{Int,1}
				VkSS     	::Array{Int,1}
				$(fieldsSS_controls...)
				distr_m		::Array{Int,1}
				distr_k		::Array{Int,1}
				distr_y		::Array{Int,1}
				D     	::Array{Int,1}
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
	Copula
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
  