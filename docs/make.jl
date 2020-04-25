push!(LOAD_PATH,"../src/")
using Documenter, HANKEstim
makedocs(sitename="Documentation for HANKEstim module",
            pages=[
                "Home" => "index.md",
                "Steady state" => "steadystate.md",
                "Linearization" => "linearization.md",
                "Estimation" => "estimation.md"
            ], format = Documenter.HTML(prettyurls = false))
