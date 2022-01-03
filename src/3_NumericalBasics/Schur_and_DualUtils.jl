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


function tot_dual(x::ForwardDiff.Dual)
    a = sum(ForwardDiff.partials(x,:))
    return a
end
function realpart(x::ForwardDiff.Dual)
    a = ForwardDiff.value(x)
    return a
end
function realpart(x::Float64)
    a = x
    return a
end
function dualpart(x::ForwardDiff.Dual)
    a = ForwardDiff.partials.(x)
    b = a.values[1]

    return b
end
function dualpart(x::Float64)

    b = x

    return b
end