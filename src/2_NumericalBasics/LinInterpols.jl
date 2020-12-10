#-------------------------------------------------------------------------------
## Linear interpolations ##
#-------------------------------------------------------------------------------
# of one function
function mylinearinterpolate(xgrd::AbstractVector,#äUnion{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd::AbstractVector,#äUnion{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    xeval::AbstractVector)#äUnion{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}})
    ind = zeros(Int,length(xeval))
    n_xgrd = length(xgrd)
    yeval= Array{eltype(ygrd),1}(undef,length(xeval))

    n_xgrd = length(xgrd)
    @views for i in eachindex(xeval)
        xi = xeval[i]
        if xi .> xgrd[n_xgrd-1]
            iL = n_xgrd - 1
        elseif xi .< xgrd[2]
            iL = 1
        else
            iL      = locate(xi, xgrd)
        end
        iR = iL+1
        xL = xgrd[iL]
        wR = (xi .- xL)./ (xgrd[iR] .- xL)
        wL = 1.0-wR
        yeval[i] = wL .* ygrd[iL] .+ wR .* ygrd[iR]
    end

    return yeval
end

# of two functions on the same grids
function mylinearinterpolate_mult2(xgrd::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd1::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd2::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    xeval::AbstractVector)#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}})
    yeval1= Array{eltype(ygrd1),1}(undef,length(xeval))
    yeval2= Array{eltype(ygrd1),1}(undef,length(xeval))

    n_xgrd = length(xgrd)
    @views @inbounds begin
        for i in eachindex(xeval)
        xi = xeval[i]
        if xi .> xgrd[n_xgrd-1]
            iL = n_xgrd - 1
        elseif xi .< xgrd[2]
            iL = 1
        else
            iL      = locate(xi, xgrd)
        end
        iR = iL+1
        xL = xgrd[iL]
        wR = (xeval[i] .- xL)./ (xgrd[iR] .- xL)
        #wL = 1.0-wR
        y1L = ygrd1[iL]
        y2L = ygrd2[iL]
        yeval1[i] = y1L .+ wR .* (ygrd1[iR] - y1L)
        yeval2[i] = y2L .+ wR .* (ygrd2[iR] - y2L)
    end
end
    return yeval1, yeval2
end

# of three functions on the same grids
function mylinearinterpolate_mult3(xgrd::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd1::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd2::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    ygrd3::AbstractVector,#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}},
    xeval::AbstractVector)#Union{AbstractArray{Float64,1},AbstractArray{DualNumbers.Dual{Float64},1}})

    yeval1= Array{eltype(ygrd1),1}(undef,length(xeval))
    yeval2= Array{eltype(ygrd2),1}(undef,length(xeval))
    yeval3= Array{eltype(ygrd3),1}(undef,length(xeval))

    n_xgrd = length(xgrd)
    for i in eachindex(xeval)
        xi = xeval[i]
        if xi .> xgrd[n_xgrd-1]
            iL = n_xgrd - 1
        elseif xi .< xgrd[2]
            iL = 1
        else
            iL      = locate(xi, xgrd)
        end
        iR = iL+1
        xL = xgrd[iL]
        wR = (xi .- xL)./ (xgrd[iR] .- xL)
        y1L = ygrd1[iL]
        y2L = ygrd2[iL]
        y3L = ygrd3[iL]
        yeval1[i] = y1L .+ wR .* (ygrd1[iR] - y1L)
        yeval2[i] = y2L .+ wR .* (ygrd2[iR] - y2L)
        yeval3[i] = y3L .+ wR .* (ygrd3[iR] - y3L)
    end
    return yeval1, yeval2, yeval3
end

#---------------- Billinear interpolation -------------------------------------
function mylinearinterpolate2(xgrd1::AbstractVector, xgrd2::AbstractVector, ygrd::AbstractArray,
    xeval1::AbstractVector, xeval2::AbstractVector)

    yeval   = zeros(eltype(xeval1),length(xeval1),length(xeval2))
    n_xgrd1 = length(xgrd1)
    n_xgrd2 = length(xgrd2)
    weight1 = Array{eltype(xeval1),1}(undef,length(xeval1))
    weight2 = Array{eltype(xeval2),1}(undef,length(xeval2))
    ind1    = zeros(Int,length(xeval1))
    ind2    = zeros(Int,length(xeval2))
    @views for i in eachindex(xeval1)
        xi = xeval1[i]
        if xi .> xgrd1[n_xgrd1-1]
            iL = n_xgrd1 - 1
        elseif xi .< xgrd1[2]
            iL = 1
        else
            iL      = locate(xi, xgrd1)
        end
        ind1[i]      = copy(iL)
        weight1[i]   = copy((xi .- xgrd1[iL])./ (xgrd1[iL.+1] .- xgrd1[iL]))
    end

    @views for i in eachindex(xeval2)
        xi = xeval2[i]
        if xi .> xgrd2[n_xgrd2-1]
            iL = n_xgrd2 - 1
        elseif xi .< xgrd2[2]
            iL = 1
        else
            iL      = locate(xi, xgrd2)
        end
        ind2[i]      = copy(iL)
        weight2[i]   = copy((xi .- xgrd2[iL])./ (xgrd2[iL.+1] .- xgrd2[iL]))
    end

    for j in eachindex(xeval2)
        w2R = weight2[j]
        w2L = 1.0-w2R
        for i in eachindex(xeval1)
            w1R = weight1[i]
            w1L = 1.0-w1R
            aux = w2L*(w1L*ygrd[ind1[i],ind2[j]] + w1R*ygrd[ind1[i]+1,ind2[j]]) +
                  w2R*(w1L*ygrd[ind1[i],ind2[j]+1] + w1R*ygrd[ind1[i]+1,ind2[j]+1])
            yeval[i,j] = aux[1]
        end
    end

    return yeval
end

@doc raw"""
    mylinearinterpolate3(xgrd1,xgrd2,xgrd3,ygrd,xeval1,xeval2,xeval3)

Trilineary project `ygrd` on (`xgrd1`,`xgrd2`,`xgrd3`) and use it to
interpolate value at (`xeval1`,`xeval2`,`xeval3`).

# Example
```jldoctest
julia> xgrd = [1.0,6.0];
julia> f((x,y,z)) = x+y+z;
julia> ygrd = f.(collect(Iterators.product(xgrid,xgrid,xgrid));
julia> xeval = [3.0,5.0];
julia> mylinearinterpolate3(xgrd,xgrd,xgrd,ygrd,xeval,xeval,xeval)
2x2x2 Array{Float64,3}:
[:,:,1] =
 9.0 11.0
11.0 13.0
[:,:,2] =
11.0 13.0
13.0 15.0
```
"""
function mylinearinterpolate3(xgrd1::AbstractVector, xgrd2::AbstractVector,
    xgrd3::AbstractVector, ygrd::AbstractArray,
    xeval1::AbstractVector, xeval2::AbstractVector, xeval3::AbstractVector)
    ## Define variables ##
    yeval   = Array{eltype(xeval1),3}(undef,(length(xeval1),length(xeval2),length(xeval3)))
    n_xgrd1 = length(xgrd1)
    n_xgrd2 = length(xgrd2)
    n_xgrd3 = length(xgrd3)
    weight1r = Array{eltype(xeval1),1}(undef,length(xeval1))
    weight2r = Array{eltype(xeval2),1}(undef,length(xeval2))
    weight3r = Array{eltype(xeval3),1}(undef,length(xeval3))
    ind1    = Array{Int,1}(undef,length(xeval1))
    ind2    = Array{Int,1}(undef,length(xeval2))
    ind3    = Array{Int,1}(undef,length(xeval3))
    ## find left element of xeval1 on xgrd1
    for i in eachindex(xeval1)
        xi = xeval1[i]
        if xi .>= xgrd1[n_xgrd1-1]
            iL = n_xgrd1 - 1
        elseif xi .< xgrd1[2]
            iL = 1
        else
            iL      = locate(xi, xgrd1)
        end
        ind1[i]      = copy(iL)
        xL          = xgrd1[iL]
        weight1r[i]  = copy((xi .- xL)./ (xgrd1[iL.+1] .- xL))
    end

    ## find left element of xeval2 on xgrd2
    for i in eachindex(xeval2)
        xi = xeval2[i]
        if xi .>= xgrd2[n_xgrd2-1]
            iL = n_xgrd2 - 1
        elseif xi .< xgrd2[2]
            iL = 1
        else
            iL      = locate(xi, xgrd2)
        end
        ind2[i]      = copy(iL)
        xL           = xgrd2[iL]
        weight2r[i]  = copy((xi .- xL)./ (xgrd2[iL.+1] .- xL))
    end

    ## find left element of xeval3 on xgrd3
    for i in eachindex(xeval3)
        xi = xeval3[i]
        if xi .>= xgrd3[n_xgrd3-1]
            iL = n_xgrd3 - 1
        elseif xi .< xgrd3[2]
            iL = 1
        else
            iL      = locate(xi, xgrd3)
        end
        ind3[i]      = copy(iL)
        xL           = xgrd3[iL]
        weight3r[i]  = copy((xi .- xL)./ (xgrd3[iL.+1] .- xL))
    end
@inbounds  begin
    for l in eachindex(xeval3)
        i3  = ind3[l]
        w3r = weight3r[l]
        w3l = 1.0 - w3r
        for j in eachindex(xeval2)
            i2 =ind2[j]
            w2r = weight2r[j]
            w2l = 1.0 - w2r
            for k in eachindex(xeval1)
                i1  = ind1[k]
                w1r = weight1r[k]
                w1l = 1.0 - w1r
               yeval[k,j,l] =  w1l.*(w2l.*(w3l.*ygrd[i1,i2,i3]     + w3r.*ygrd[i1,i2,i3+1])    +
                                      w2r.*(w3l.*ygrd[i1,i2+1,i3]   + w3r.*ygrd[i1,i2+1,i3+1])) +
                                w1r.*(w2l.*(w3l.*ygrd[i1+1,i2,i3]   + w3r.*ygrd[i1+1,i2,i3+1])  +
                                      w2r.*(w3l.*ygrd[i1+1,i2+1,i3] + w3r.*ygrd[i1+1,i2+1,i3+1]))
            end
        end
    end
end
    return yeval
end
