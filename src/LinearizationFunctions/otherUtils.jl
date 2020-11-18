##########################################################
# Matrix to remove one degree of freedom from distribution
#---------------------------------------------------------
function shuffleMatrix(distr, n_par)
    distr_m = sum(sum(distr,dims=3),dims=2)./sum(distr[:])
    distr_k = sum(sum(distr,dims=3),dims=1)./sum(distr[:])
    distr_y = sum(sum(distr,dims=2),dims=1)./sum(distr[:])
    Γ    = Array{Array{Float64,2},1}(undef,3)
    Γ[1] = zeros(Float64,(n_par.nm,n_par.nm-1))
    Γ[2] = zeros(Float64,(n_par.nk,n_par.nk-1))
    Γ[3] = zeros(Float64,(n_par.ny,n_par.ny-1))
    for j=1:n_par.nm-1
        Γ[1][:,j] = -distr_m[:]
        Γ[1][j,j] = 1-distr_m[j]
        Γ[1][j,j] = Γ[1][j,j] - sum(Γ[1][:,j])
    end
    for j=1:n_par.nk-1
        Γ[2][:,j] = -distr_k[:]
        Γ[2][j,j] = 1-distr_k[j]
        Γ[2][j,j] = Γ[2][j,j] - sum(Γ[2][:,j])
    end
    for j=1:n_par.ny-1
        Γ[3][:,j] = -distr_y[:]
        Γ[3][j,j] = 1-distr_y[j]
        Γ[3][j,j] = Γ[3][j,j] - sum(Γ[3][:,j])
    end

    return Γ
end
##########################################################
# Schur Decomposition
#---------------------------------------------------------
function complex_schur(A, B)
    F = LinearAlgebra.schur(complex(A), complex(B))
    α::Vector{complex(promote_type(eltype(A), eltype(B)))} = F.alpha
    λ = abs.(α) ./ abs.(F.beta)
    select_ev = λ .>= 1.0
    # select_ev = abs.(λ) .>= 1.0
    nk  = sum(select_ev) # Number of state Variables based on Eigenvalues
    return F, select_ev, nk, λ
end
##########################################################
# Indexes
#---------------------------------------------------------

@doc raw"""
	@make_fnaggr(fn_name,aggr_names)

Create function `fn_name` that returns an instance of `struct` `IndexStructAggr`
(created by [`@make_struct_aggr`](@ref)), mapping aggregate states and controls to values
`1` to `length(aggr_names)` (both steady state and deviation from it).
"""
macro make_fnaggr(fn_name,  aggr_names)
	state_names=Symbol.(eval((aggr_names)))
	n_states = length(state_names)

	fieldsSS_states = [:( $i) for i = 1:n_states]
	fields_states = [:($i) for i = 1:n_states]
	esc(quote
		function $(fn_name)(n_par)
		    indexes = IndexStructAggr(
				$(fieldsSS_states...),
				$(fields_states...)
				)
			return indexes
		end
	end)
end

@doc raw"""
	@make_fn(fn_name,s_names,c_names)

Create function `fn_name` that returns an instance of `struct` `IndexStruct`
(created by [`@make_struct`](@ref)), mapping states and controls to indexes
inferred from numerical parameters and compression indexes.
"""
macro make_fn(fn_name,  s_names, c_names)
    # fields=[:($(entry.args[1])::$(entry.args[2])) for entry in var_names]
	# fieldsSS=[:($(Symbol((entry.args[1]), "SS"))::$(entry.args[2])) for entry in var_names]
	state_names=Symbol.(eval((s_names)))
	n_states = length(state_names)
	control_names=Symbol.(eval((c_names)))
	n_controls = length(control_names)

	fieldsSS_states = [:((tNo + tNo2) + $i) for i = 1:n_states]
	fields_states = [:(tNo + tNo4 - 3 + $i) for i = 1:n_states]
	fieldsSS_controls = [:(tNo + 3 * tNo2 + $i) for i = n_states .+ (1:n_controls)]
	fields_controls = [:(tNo + tNo3 + tNo4 - 3 + $i) for i = n_states .+ (1:n_controls)]
	esc(quote
		function $(fn_name)(n_par, compressionIndexesVm, compressionIndexesVk, compressionIndexesD)
		    tNo = n_par.nm + n_par.nk + n_par.ny
		    tNo2 = n_par.nm * n_par.nk * n_par.ny
		    tNo3 = length(compressionIndexesVm) + length(compressionIndexesVk)
			tNo4 = length(compressionIndexesD)
		    indexes = IndexStruct(
		        1:n_par.nm, # distr_m_SS
		        (n_par.nm+1):(n_par.nm+n_par.nk), # distr_k_SS
		        (n_par.nm+n_par.nk+1):(tNo), # distr_y_SS
				(tNo + 1):(tNo + tNo2), # VDSS
				$(fieldsSS_states...),
				((tNo + tNo2) + $(n_states) + 1):((tNo + tNo2) + tNo2 + $(n_states)), # VmSS
		        ((tNo + tNo2) + tNo2 + $(n_states) + 1):((tNo + tNo2) + 2 * tNo2 + $(n_states)), # VkSS
				$(fieldsSS_controls...),
				1:(n_par.nm - 1), # distr_m
		        n_par.nm:(n_par.nm + n_par.nk - 2), # distr_k
		        (n_par.nm + n_par.nk - 1):(tNo - 3), # distr_y
				(tNo - 2):(tNo + tNo4 - 3), # VD
				$(fields_states...),
				(tNo + tNo4 + $(n_states) - 2):(tNo + tNo4 + length(compressionIndexesVm) + $(n_states) - 3), # Vm
		        (tNo + tNo4 + length(compressionIndexesVm) + $(n_states) -2):(tNo + tNo4 + length(compressionIndexesVm) + length(compressionIndexesVk) + $(n_states) - 3), # Vk
				$(fields_controls...)
				)
			return indexes
		end
	end)
end


@doc raw"""
	@make_deriv(n_FD_s)

Set dual number seeds at chunk of size `n_FD_s` (starting at index `i`)
of the argument `X` (`XPrime`) in [`Fsys()`](@ref) and save matrix
of partial derivatives (for this group of variables) in `Deriv(i)` (`DerivPrime(i)`).

Require `length_X0` and remaining arguments of [`Fsys()`](@ref) in namespace.
"""
macro make_deriv(n_FD_s)
	n_FD = eval(n_FD_s)
	a = Matrix{Float64}(I, eval(n_FD), eval(n_FD))
	fn_rows =Expr(:call, :.+, [:(((1:length_X0) .== i + ($j-1)) * ForwardDiff.Dual(0.0, tuple($(a[j, :])...))) for j = 1:n_FD]...)
	esc(quote
    	Deriv(i)      = ForwardDiff.partials.(Fsys(X0 .+ $(fn_rows), X0, XSS, m_par,
                                n_par, indexes, Γ, compressionIndexes,DC,IDC, DCD, IDCD),:)
    	DerivPrime(i) = ForwardDiff.partials.(Fsys(X0, X0 .+ $(fn_rows), XSS, m_par,
                                n_par, indexes, Γ, compressionIndexes,DC,IDC, DCD, IDCD),:)
		end
	)
end

@doc raw"""
	@make_deriv_estim(n_FD_s)

Set dual number seeds at chunk of size `n_FD_s` (starting at index `i`)
of the argument `X` (`XPrime`) in [`Fsys_agg()`](@ref) and save matrix
of partial derivatives (for this group of variables) in `Deriv(i)` (`DerivPrime(i)`).

Require `length_X0` and remaining arguments of [`Fsys_agg()`](@ref) in namespace.
"""
macro make_deriv_estim(n_FD_s)
	n_FD = eval(n_FD_s)
	a = Matrix{Float64}(I, eval(n_FD), eval(n_FD))
	fn_rows =Expr(:call, :.+, [:(((1:length_X0) .== i + ($j-1)) * ForwardDiff.Dual(0.0, tuple($(a[j, :])...))) for j = 1:n_FD]...)
	esc(quote
    	Deriv(i)      = ForwardDiff.partials.(Fsys_agg(X0 .+ $(fn_rows), X0, XSSaggr, distrSS, m_par,
                                n_par, indexes_aggr),:)
    	DerivPrime(i) = ForwardDiff.partials.(Fsys_agg(X0, X0 .+ $(fn_rows), XSSaggr, distrSS, m_par,
                                n_par, indexes_aggr),:)
		end
	)
end
