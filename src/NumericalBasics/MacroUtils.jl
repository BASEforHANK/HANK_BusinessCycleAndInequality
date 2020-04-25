@doc raw"""
    @generate_equations(var_names)

Write out the expansions around steady state for all variables in `var_names`,
i.e. generate code that reads aggregate states/controls from steady state deviations.

Equations take the form of (with variable `r` as example):
- `r       = exp.(Xss[indexes.rSS] .+ X[indexes.r])`
- `rPrime  = exp.(Xss[indexes.rSS] .+ XPrime[indexes.r])`
"""
macro generate_equations(var_names)
    ex = quote end # initialize expression to append to
    var_names=eval(var_names)
    for j in var_names # loop over variables to generate
        i = Symbol(j)
        varnamePrime = Symbol(i,"Prime")
        varnameSS = Symbol(i,"SS")
        ex_aux = quote
            $i = exp.(Xss[indexes.$varnameSS] .+ X[indexes.$i])
            $varnamePrime = exp.(Xss[indexes.$varnameSS] .+ XPrime[indexes.$i])
        end

        append!(ex.args, ex_aux.args) # append to expression
    end

    return esc(ex)
end

"""
    @writeXSS()

Write all single steady state variables into vectors XSS / XSSaggr.

# Requires
globals `state_names`, `control_names`, `aggr_names`
"""
macro writeXSS()
        ex = quote
                XSS =   [ distr_m_SS[:];distr_k_SS[:]; distr_y_SS[:]]
        end
        for j in state_names
                varnameSS = Symbol(j,"SS")
                ex_aux = quote
                    append!(XSS,log($varnameSS))
                end
                append!(ex.args, ex_aux.args)
        end

        ex_aux= quote
                append!(XSS,[VmSS[:]; VkSS[:]]) # value function controls
        end
        append!(ex.args, ex_aux.args)
        for j in control_names
            varnameSS = Symbol(j,"SS")
            ex_aux = quote
                append!(XSS,log($varnameSS))
            end
            append!(ex.args, ex_aux.args)
        end
        ex_aux= quote XSSaggr=[0.0] end
        append!(ex.args, ex_aux.args)
        for j in aggr_names
                varnameSS = Symbol(j,"SS")
                ex_aux = quote
                    append!(XSSaggr,log($varnameSS))
                end
                append!(ex.args, ex_aux.args)
        end
        ex_aux= quote deleteat!(XSSaggr,1) end
        append!(ex.args, ex_aux.args)

    return esc(ex)
end

macro include(filename::AbstractString)
    path = joinpath(dirname(String(__source__.file)), filename)
    return esc(Meta.parse("quote; " * read(path, String) * "; end").args[1])
end
