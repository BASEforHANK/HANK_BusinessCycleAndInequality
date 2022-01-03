@doc raw"""
    EGM_policyupdate(EVm,EVk,Qminus,πminus,RBminus,Tshock,inc,n_par,m_par,warnme)

Find optimal policies, given marginal continuation values `EVm`, `EVk`, today's
prices [`Qminus`, `πminus`,`RBminus`], and income [`inc`], using the
Endogenous Grid Method.

Optimal policies are defined on the fixed grid, but optimal asset choices (`m` and `k`)
are off-grid values.

# Returns
- `c_a_star`,`m_a_star`,`k_a_star`,`c_n_star`,`m_n_star`: optimal (on-grid) policies for
    consumption [`c`], liquid [`m`] and illiquid [`k`] asset, with [`a`] or
    without [`n`] adjustment of illiquid asset
"""
function EGM_policyupdate(EVm::Array,
                          EVk::Array,
                          Qminus::Real,
                          πminus::Real,
                          RBminus::Real,
                          Tshock::Real,
                          inc::Array,
                          n_par::NumericalParameters,
                          m_par::ModelParameters, 
                          warnme::Bool)

################### Copy/read-out stuff#####################################
β::Float64  = m_par.β
borrwedge   = m_par.Rbar.*Tshock
# inc[1] = labor income , inc[2] = rental income,
# inc[3]= liquid assets income, inc[4] = capital liquidation income
inc_lab     = inc[1]
inc_rent    = inc[2]
inc_LA      = inc[3]
inc_IA      = inc[4]
n           = size(EVm)
mmax        = n_par.grid_m[end]
kmax        = n_par.grid_k[end]

############################################################################
## EGM Step 1: Find optimal liquid asset holdings in the constrained case ##
############################################################################
EMU         = EVm .* β
c_star_n    = invmutil(EMU) # 6% of time with rolled out power function

# Calculate assets consistent with choices being [m']
# Calculate initial money position from the budget constraint
# that leads to the optimal consumption choice
m_star_n    = (c_star_n .+ n_par.mesh_m .- inc_lab .- inc_rent)
# Apply correct interest rate
m_star_n   .= m_star_n ./ (RBminus ./ πminus .+ borrwedge.*(m_star_n.<0))  # apply borrowing rate

# Next step: Interpolate w_guess and c_guess from new k-grids
# using c[s,h,m"], m(s,h,m")
# Interpolate grid().m and c_n_aux defined on m_star_n over grid().m

# Check monotonicity of m_star_n
if warnme
    m_star_aux    = reshape(m_star_n,(n[1], n[2]*n[3]))
    if any(any(diff(m_star_aux, dims=1).<0))
        @warn "non monotone future liquid asset choice encountered"
    end
end

# Policies for tuples (c*,m*,y) are now given. Need to interpolate to return to
# fixed grid.
c_n_star    = Array{eltype(c_star_n),3}(undef,(n[1],n[2], n[3]))#zeros(eltype(c_star_n), size(c_star_n)) # Initialize c_n-container
m_n_star    = Array{eltype(c_star_n),3}(undef,(n[1],n[2], n[3]))#zeros(eltype(c_star_n), size(c_star_n)) # Initialize m_n-container

@inbounds @views  begin
    for jj=1:n[3] # Loop over income states
        for kk = 1:n[2] # Loop over capital states
            cc, mn = mylinearinterpolate_mult2(m_star_n[:,kk,jj], c_star_n[:,kk,jj],n_par.grid_m, n_par.grid_m)
            c_n_star[:,kk,jj]=cc
            m_n_star[:,kk,jj]=mn
            # Check for binding borrowing constraints, no extrapolation from grid
            bcpol = m_star_n[1,kk,jj]
            for mm= 1:n[1]
                if n_par.mesh_m[mm,kk,jj] .<bcpol
                    c_n_star[mm,kk,jj] = inc_lab[mm,kk,jj] .+ inc_rent[mm,kk,jj] .+ inc_LA[mm,kk,jj] .- n_par.grid_m[1]
                    m_n_star[mm,kk,jj] = n_par.grid_m[1]
                end
                if mmax  .< m_n_star[mm,kk,jj]
                    m_n_star[mm,kk,jj] = mmax
                end
            end
        end
    end
end

#-------------------------END OF STEP 1-----------------------------

############################################################################
## EGM Step 2: Find Optimal Portfolio Combinations                        ##
############################################################################


# Find an m_a* for given k' that yield the same expected future marginal value
# for liquid and illiquid assets:
term1           = (β/Qminus) * EVk                      # expected marginal value of illiquid investment
E_return_diff   = term1 .- EMU                          # difference conditional on future asset holdings on grid
m_a_aux1        = Fastroot(n_par.grid_m,E_return_diff)  # Find indifferent m by interpolation of two neighboring points a, b ∈ grid_m with:  E_return_diff(a) < 0 < E_return_diff(b)
                                                        # (Fastroot does not allow for extrapolation and uses non-negativity constraint and monotonicity)
m_a_aux         = reshape(m_a_aux1, (n[2],n[3]))

###########################################################################
## EGM Step 3: Constraints for money and capital are not binding         ##
###########################################################################
# Interpolation of psi()-function at m*_n[m,k]
aux_index       = (0:(n[2]*n[3])-1)*n[1]                        # auxiliary to move to linear indexing
EMU_star        = Array{eltype(m_a_aux),2}(undef,(n[2],n[3]))   # container
step            = diff(n_par.grid_m)                            # Stepsize on grid()

# Interpolate EMU[m",k',s'*h',M',K'] over m*_n[k"], m-dim is dropped
for j in eachindex(m_a_aux)
    xi          = m_a_aux[j]
    # find indexes on grid next smallest to optimal policy
    if xi.> n_par.grid_m[n[1]-1]                                # policy is larger than highest grid point 
        idx     = n[1]-1
    elseif xi.<= n_par.grid_m[1]                                # policy is smaller than lowest grid point
        idx     = 1
    else
        idx     = locate(xi,n_par.grid_m)                       # use exponential search to find grid point closest to policy (next smallest)
    end
    
    s           = (xi .- n_par.grid_m[idx])./step[idx]          # Distance of optimal policy to next grid point

    EMU_star[j] = EMU[idx .+ aux_index[j]].*(1.0 -s) .+
                        s.*(EMU[idx .+ aux_index[j].+1])        # linear interpolation
end


c_a_aux         = invmutil(EMU_star)

# Resources that lead to capital choice
# k'= c + m*(k") + k" - w*h*N
# = value of todays cap and money holdings
Resource        = c_a_aux .+ m_a_aux .+ inc_IA[1,:,:] .- inc_lab[1,:,:]

# Money constraint is not binding, but capital constraint is binding
m_star_zero     = m_a_aux[1,:] # Money holdings that correspond to k'=0:  m*(k=0)

# Use consumption at k"=0 from constrained problem, when m" is on grid()
aux_c           = reshape(c_star_n[:,1,:],(n[1], n[3]))
aux_inc         = reshape(inc_lab[1,1,:],(1, n[3]))
cons_list       = Array{Array{eltype(c_star_n)}}(undef,n[3],1)
res_list        = Array{Array{eltype(c_star_n)}}(undef,n[3],1)
mon_list        = Array{Array{eltype(c_star_n)}}(undef,n[3],1)
cap_list        = Array{Array{eltype(c_star_n)}}(undef,n[3],1)

for j=1:n[3]
    # When choosing zero capital holdings, HHs might still want to choose money
    # holdings smaller than m*(k'=0)
    if m_star_zero[j]>n_par.grid_m[1]
        # Calculate consumption policies, when HHs chooses money holdings
        # lower than m*(k"=0) and capital holdings k"=0 and save them in cons_list
        log_index    = n_par.grid_m.<m_star_zero[j]
        # aux_c is the consumption policy under no cap. adj. (fix k=0), for m<m_a*(k'=0)
        c_k_cons     = aux_c[log_index,j]
        cons_list[j] = c_k_cons; #Consumption at k"=0, m"<m_a*(0)
        # Required Resources: Money choice + Consumption - labor income
        # Resources that lead to k"=0 and m'<m*(k"=0)
        res_list[j]  = n_par.grid_m[log_index] .+ c_k_cons .- aux_inc[j]
        mon_list[j]  = n_par.grid_m[log_index]
        cap_list[j]  = zeros(eltype(EVm),sum(log_index))
    else
        cons_list[j] = zeros(eltype(EVm),0) #Consumption at k"=0, m"<m_a*(0)
        # Required Resources: Money choice + Consumption - labor income
        # Resources that lead to k"=0 and m'<m*(k"=0)
        res_list[j]  = zeros(eltype(EVm),0)
        mon_list[j]  = zeros(eltype(EVm),0)
        cap_list[j]  = zeros(eltype(EVm),0)
    end
end

# Merge lists
c_a_aux             = reshape(c_a_aux,(n[2], n[3]))
m_a_aux             = reshape(m_a_aux,(n[2], n[3]))

for j=1:n[3]
    cons_list[j]    = append!(cons_list[j], c_a_aux[:,j])
    res_list[j]     = append!(res_list[j],  Resource[:,j])
    mon_list[j]     = append!(mon_list[j],  m_a_aux[:,j])
    cap_list[j]     = append!(cap_list[j],  n_par.grid_k)
end

####################################################################
## EGM Step 4: Interpolate back to fixed grid                     ##
####################################################################
c_a_star            = Array{eltype(c_star_n),3}(undef,(n[1],n[2], n[3]))
m_a_star            = Array{eltype(c_star_n),3}(undef,(n[1],n[2], n[3]))#zeros(eltype(c_star_n), (n[1],n[2], n[3]))
k_a_star            = Array{eltype(c_star_n),3}(undef,(n[1],n[2], n[3]))#zeros(eltype(c_star_n), (n[1],n[2], n[3]))
Resource_grid       = reshape(inc_IA .+ inc_LA .+ inc_rent,(n[1].*n[2], n[3]))
labor_inc_grid      = inc_lab[1,1,:][:]#reshape(inc_lab,(n[1]*n[2], n[3]))
log_index2          = zeros(Bool,n[1].*n[2])

@views @inbounds begin
    for j=1:n[3]
        # Check monotonicity of resources
        if warnme
            if any(diff(res_list[j]).<0)
                @warn "non monotone resource list encountered"
            end
        end
        # when at most one constraint binds:
        log_index2[:] = reshape(Resource_grid[:,j],n[1]*n[2]).<res_list[j][1]
        # Lowest value of res_list corresponds to m_a"=0 and k_a"=0.
        c_a_star1, m_a_star1, k_a_star1 = mylinearinterpolate_mult3(res_list[j],cons_list[j],mon_list[j],cap_list[j],Resource_grid[:,j])

        # Any resources on grid smaller then res_list imply that HHs consume all
        # resources plus income.
        # When both constraints are binding:
         c_a_star1[log_index2]  = Resource_grid[log_index2,j] .+ labor_inc_grid[j] .- n_par.grid_m[1]
         m_a_star1[log_index2] .= n_par.grid_m[1]
         k_a_star1[log_index2] .= 0.0
         # Check if policies go beyond largest gr
         for kk = 1:n[2]               
            for mm = 1:n[1]
                runind                  = mm + (kk-1)*n[1]  # linear indexing
                mp                      = m_a_star1[runind] # liquid asset saving
                kp                      = k_a_star1[runind] # illiquid asset saving  
                c_a_star[mm,kk,j]       = c_a_star1[runind] # consumption policy
                if mp <mmax
                    m_a_star[mm,kk,j]   = mp
                else
                    m_a_star[mm,kk,j]   = mmax
                end
                if kp < kmax
                    k_a_star[mm,kk,j]   = kp
                else
                    k_a_star[mm,kk,j]   = kmax
                end
            end
        end
    end
end

    return c_a_star, m_a_star, k_a_star, c_n_star, m_n_star
end

