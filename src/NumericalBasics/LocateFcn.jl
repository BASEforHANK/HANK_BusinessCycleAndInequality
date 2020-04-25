#-------------------------------------------------------------------------------
## Locate function ##
#-------------------------------------------------------------------------------
locate(x::Number, xx::AbstractVector) = exp_search(x, xx)
function bin_search(x::Number, xx::AbstractVector)
    # real(8), intent(in), dimension(N)		:: xx  ! lookup table
    # integer, intent(in)						:: N   ! no elements of lookup table
    # real(8), intent(in)						:: x   ! Value whose nearest neighbors are to be found
    # integer, intent(out)					:: j   ! returns j if xx(j)<x<xx(j+1),
    #												   ! 0 if value is to the left of grid,
    #												   ! N if value is to the right of grid

    N = length(xx)
    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        jl = 1
        ju = N
        while (ju-jl) != 1
            jm = div((ju+jl),2)

@inbounds   if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end

    return j
end

function exp_search(x::Number, xx::AbstractVector)
    # real(8), intent(in), dimension(N)		:: xx  ! lookup table
    # integer, intent(in)						:: N   ! no elements of lookup table
    # real(8), intent(in)						:: x   ! Value whose nearest neighbors are to be found
    # integer, intent(out)					:: j   ! returns j if xx(j)<x<xx(j+1),
    #												   ! 0 if value is to the left of grid,
    #												   ! N if value is to the right of grid
@inbounds begin
    N = length(xx)
    if x <= xx[1]
        j = 1
    elseif x >= xx[N]
        j = N
    else
        bound = 1
        while bound < N && x>xx[bound]
            bound *=2
        end
        jl = div(bound,2)
        ju = min(N,bound)
        while (ju-jl) != 1
            jm = div((ju+jl),2)

            if x .> xx[jm]
                jl = jm
            else
                ju = jm
            end
        end
        j = jl
    end
end
    return j
end
