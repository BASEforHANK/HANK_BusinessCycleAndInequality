#---------------------------------------------------------
# Integration (Gauss Chebychev)
#---------------------------------------------------------
function my_integrate(f::Function, a::Number, b::Number)
    nodes, weights = my_qnwcheb(500, a, b)
    I = weights'*f.(nodes)
end

function my_qnwcheb(n::Integer, a::Number, b::Number)
    nodes = (b + a) / 2 .- (b - a) / 2 .* cos.(pi / n .* (0.5:(n - 0.5)))
    weights = ((b - a) / n) .* (cos.(pi / n .* ((1:n) .- 0.5) * (2:2:n-1)') *
    (-2.0 ./ ((1:2:n-2) .* (3:2:n))) .+ 1)
    return nodes, weights
end
