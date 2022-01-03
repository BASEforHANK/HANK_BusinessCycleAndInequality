function cdf_to_pdf(cdf)
    Y = copy(cdf)
    for j=1:length(size(cdf))
        t = collect(size(cdf))
        t[j] = 1
        aux = zeros(eltype(Y),Tuple(t))
        Y = diff(cat([aux, Y]...;dims=j); dims = j)
    end
    return Y
end

function pdf_to_cdf(pdf)
    Y = copy(pdf)
    for j=1:length(size(pdf))
        Y = cumsum(Y; dims = j)
    end
    return Y
end