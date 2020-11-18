#----------------------------------------------------------------------------
# Basic Functions: Return on capital, marginal utility and its inverse
#---------------------------------------------------------------------------
mutil(c::Array, ξ::Float64)      = 1.0 ./ ((c.*c).*(c.*c)) # c.^ξ#
invmutil(mu::Array, ξ::Float64)  = 1.0 ./ (sqrt.(sqrt.(mu))) # mu.^(1.0./ξ)#

# Incomes (K:capital, A: TFP): Interest rate = MPK.-δ, Wage = MPL, profits = Y-wL-(r+\delta)*K
interest(K::Number, A::Number,N::Number, m_par::ModelParameters) = A .* m_par.α .* (K ./ N) .^(m_par.α - 1.0) .- m_par.δ_0
wage(K::Number, A::Number,N::Number, m_par::ModelParameters)     = A .* (1 - m_par.α) .* (K./N) .^m_par.α
output(K::Number, A::Number,N::Number, m_par::ModelParameters)   = A .* K .^(m_par.α) .* N .^(1 - m_par.α)
employment(K::Number, A::Number, m_par::ModelParameters)         = (A .* (1.0 - m_par.α) .* (m_par.τ_lev .* (1.0 - m_par.τ_prog)).^(1.0 / (1.0 - m_par.τ_prog)) 
                                                                    .* K .^(m_par.α )).^((1.0 - m_par.τ_prog)./(m_par.γ + m_par.τ_prog + (m_par.α) .* (1 - m_par.τ_prog)))

@doc raw"""
    distrSummaries(distr,c_a_star,c_n_star,n_par,inc,incgross,m_par)

Compute distributional summary statistics, e.g. Gini indexes, top-10%
income and wealth shares, and 10%, 50%, and 90%-consumption quantiles.

# Arguments
- `distr`: joint distribution over bonds, capital and income ``(m \times k \times y)``
- `c_a_star`,`c_n_star`: optimal consumption policies with [`a`] or without [`n`]
    capital adjustment
- `n_par::NumericalParameters`, `m_par::ModelParameters`
- `inc`: vector of (on grid-)incomes, consisting of labor income (scaled by ``\frac{\gamma-\tau^P}{1+\gamma}``, plus labor union-profits),
    rental income, liquid asset income, capital liquidation income,
    labor income (scaled by ``\frac{1-\tau^P}{1+\gamma}``, without labor union-profits),
    and labor income (without scaling or labor union-profits)
- `incgross`: vector of (on grid-) *pre-tax* incomes, consisting of
    labor income (without scaling, plus labor union-profits), rental income,
    liquid asset income, capital liquidation income,
    labor income (without scaling or labor union-profits)
"""
function distrSummaries(distr::AbstractArray,c_a_star::AbstractArray,
                        c_n_star::AbstractArray, n_par::NumericalParameters,
                        inc::AbstractArray,incgross::AbstractArray, m_par::ModelParameters)
    ## Distributional summaries
    mplusk = zeros(n_par.nk.*n_par.nm)
    for k = 1:n_par.nk
        for m = 1:n_par.nm
            mplusk[m+(k-1)*n_par.nm] = n_par.grid_m[m]+n_par.grid_k[k];
        end
    end
    IX=sortperm(mplusk)
    mplusk          = mplusk[IX]
    moneycapital_pdf= sum(distr, dims=3)
    moneycapital_pdf= moneycapital_pdf[IX];
    moneycapital_cdf= cumsum(moneycapital_pdf);
    S               = [0; cumsum(moneycapital_pdf.*mplusk)]
    giniwealth      = 1-(sum(moneycapital_pdf.*(S[1:end-1]+S[2:end]))/S[end]);

    distr_m = sum(distr,dims=(2,3))[:]
    distr_k = sum(distr,dims=(1,3))[:]
    distr_y = sum(distr,dims=(1,2))[:]

    share_borrower = sum(distr_m.*(n_par.grid_m.<0))

    p50             = count(moneycapital_cdf.<0.5)+1;
    p90             = count(moneycapital_cdf.<0.9)+1;
    w9050           = mplusk[p90]./mplusk[p50];
    FN_wealthshares = cumsum(mplusk.*moneycapital_pdf)./sum(mplusk.*moneycapital_pdf);
    w90share        = 1.0 - mylinearinterpolate(moneycapital_cdf,FN_wealthshares,[0.9])[1]


    x               = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    c               = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    distr_x         = zeros(eltype(c_a_star),(n_par.nm, n_par.nk,n_par.ny,2))
    x[:,:,:,1]      = c_a_star
    x[:,:,:,2]      = c_n_star;
    aux_x           = inc[5]#(1 .- m_par.τ_bar).*(1.0 ./ m_par.μw).*w.*N.*n_par.mesh_y./(m_par.γ+1);
    aux_x[:,:,end]  = zeros(n_par.nm,n_par.nk);
    c[:,:,:,1]      = x[:,:,:,1] +aux_x;
    c[:,:,:,2]      = x[:,:,:,2] +aux_x;
    distr_x[:,:,:,1] = m_par.λ .* distr
    distr_x[:,:,:,2] = (1-m_par.λ) .* distr

    IX              = sortperm(x[:]);
    x               = x[IX];
    x_pdf           = distr_x[IX];
    S               = cumsum(x_pdf.*x);
    S               = [0 S'];
    ginicompconsumption = 1-(sum(x_pdf.*(S[1:end-1]+S[2:end]))/S[end]);
    sdlogx          = sqrt(x_pdf[:]'*log.(x[:]).^2-(x_pdf[:]'*log.(x[:]))^2);

    IX              = sortperm(c[:]);
    c               = c[IX];
    c_pdf           = distr_x[IX];
    S               = cumsum(c_pdf.*c);
    c_cdf           = cumsum(c_pdf);

    p10             = count(c_cdf.<0.1)+1;
    p50             = count(c_cdf.<0.5)+1;
    p90             = count(c_cdf.<0.9)+1;

    c9010           = c[p90]./c[p10];

    p10C            = c[p10];
    p50C            = c[p50];
    p90C            = c[p90];

    S               = [0 S'];
    giniconsumption = 1-(sum(c_pdf.*(S[1:end-1]+S[2:end]))/S[end]);
    sdlogc          = sqrt(c_pdf[:]'*log.(c[:]).^2-(c_pdf[:]'*log.(c[:]))^2);

    Yidio           = inc[6]+inc[2]+inc[3] - n_par.mesh_m
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    Y_pdf           = distr[IX]
    Y_cdf           = cumsum(Y_pdf)
    p10             = count(Y_cdf.<0.1)+1
    FN_incomesharesnet = cumsum(Yidio.*Y_pdf)./sum(Yidio.*Y_pdf);
    I90sharenet        = 1.0 .- mylinearinterpolate(Y_cdf,FN_incomesharesnet,[0.9])[1]

    Yidio           = incgross[1]+incgross[2]+incgross[3] - n_par.mesh_m
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    Y_pdf           = distr[IX]
    Y_cdf           = cumsum(Y_pdf)
    FN_incomeshares = cumsum(Yidio.*Y_pdf)./sum(Yidio.*Y_pdf);
    I90share        = 1.0 .- mylinearinterpolate(Y_cdf,FN_incomeshares,[0.9])[1]

    S               = cumsum(Y_pdf.*Yidio)
    S               = [0 S']
    giniincome      = 1-(sum(Y_pdf.*(S[1:end-1]+S[2:end]))/S[end])


    Yidio           = incgross[1]
    Yidio           = Yidio[:,:,1:end-1]
    IX              = sortperm(Yidio[:])
    Yidio           = Yidio[IX]
    distr_aux       = distr[:,:,1:end-1]
    distr_aux       = distr_aux./sum(distr_aux[:])
    Y_pdf           = distr_aux[IX]
    Y_cdf           = cumsum(Y_pdf)
    p10             = count(Y_cdf.<0.1)+1
    p50             = count(Y_cdf.<0.5)+1
    y5010           = Yidio[p50]./Yidio[p10]

    sdlogy          = sqrt(Y_pdf[:]'*log.(Yidio[:]).^2-(Y_pdf[:]'*log.(Yidio[:]))^2);

    return     distr_m, distr_k, distr_y, share_borrower, giniwealth, I90share,I90sharenet, ginicompconsumption,#=
            =# sdlogx, c9010, giniconsumption, sdlogc, y5010, giniincome, sdlogy, w90share, p10C, p50C, p90C
end
