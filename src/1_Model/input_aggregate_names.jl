# This file defines the sets of aggregate shocks, states (inluding shocks), and controls
# The code checks that the number of equations in the aggregate model is equal to the number 
# of aggregate variables excluding the distributional summary statistics. The latter are not 
# contained in the aggregate model code as they are parameter free but change whenever the 
# distribution changes and do not show up in any aggregate model equation.

shock_names = [
    :A, :Z, :ZI, :μ, :μw, :Gshock, :Rshock, :Sshock, :Tprogshock
]

state_names = [
    "A", "Z", "ZI", "RB" , "μ", "μw", "σ", "Ylag",
    "Blag", "Tlag", "Ilag", "wlag", "qlag", "Clag", "av_tax_ratelag", "τproglag",
    "Gshock", "Tprogshock", "Rshock", "Sshock"
]

control_names = [
    "r", "w", "K", "π" ,"πw", "Y" ,"C", "q",  "N", "mc", "mcw", "u",
    "Ht", "av_tax_rate", "T", "I", "B","BD", "BY","TY", "mcww", "G", "τlev", "τprog", "GiniC",
    "GiniX", "TOP10Ishare", "TOP10Inetshare", "TOP10Wshare", "sdlogy",
     "Ygrowth", "Bgrowth", "Igrowth", "wgrowth", "Cgrowth", "Tgrowth", "LP", "LPXA", "unionprofits", "profits"
]
# Identify distributional summary variables (b/c no equations in aggregate model expected)
distr_names   = ["GiniC", "GiniX", "TOP10Ishare", "TOP10Inetshare", "TOP10Wshare", "sdlogy"]

# All names in one array
aggr_names          = [state_names; control_names]

# ascii names used for cases where unicode doesn't work, e.g., file saves
unicode2ascii(x)    = replace.(replace.(replace.(replace.(replace.(x,"τ"=>"tau"), "σ" => "sigma"),"π"=>"pi"),"μ"=>"mu"),"ρ"=>"rho")

state_names_ascii   = unicode2ascii(state_names)
control_names_ascii = unicode2ascii(control_names)
aggr_names_ascii    = [state_names_ascii; control_names_ascii]