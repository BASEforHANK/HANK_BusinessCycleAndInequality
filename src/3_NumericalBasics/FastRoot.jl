function Fastroot(xgrid::Vector,fx::AbstractArray)
    #fast linear interpolation root finding
    #(=one Newton step at largest negative function value)
    #   stripped down version of interp1 that accepts multiple inputs [max 3]
    #   that are interpolated over the same grids x & xi
    xgrid=xgrid[:]
    fx = reshape(fx,length(xgrid),:)
    # Make use of the fact that the difference equation is monotonically
    # increasing in m, use sign for robustness.
    n=size(fx)
    roots    = Array{eltype(fx),1}(undef,n[2])
    @views @inbounds begin
        for j=1:n[2]
            fz=view(fx,:,j)
            idx  = locate(0.0,sign.(fz))
            if idx>=n[1]
                roots[j] = xgrid[end]
            elseif idx<=1
                roots[j] = xgrid[1]
            else
                fxx  = fz[idx]
                dfxx = fz[idx.+1]  - fxx
                xx   = xgrid[idx]
                dxx  = xgrid[idx.+1] - xgrid[idx]
                roots[j] = xx-fxx.*dxx./dfxx
            end
        end
    end
    return roots
end
