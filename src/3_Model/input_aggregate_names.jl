shock_names = [
    :A, :Z, :ZI, :μ, :μw, :Gshock, :Rshock, :Sshock, :Tlevshock, :Tprogshock
]

state_names = [
    "A", "Z", "ZI", "RB" , "μ", "μw", "σ", "union_Retained", "Retained", "Ylag",
    "Blag", "Tlag", "Ilag", "wlag", "qlag", "Clag", "av_tax_ratelag", "τproglag",
    "Gshock", "Tprogshock", "Rshock", "Sshock"
]

control_names = [
    "r", "w", "K", "π" ,"πw", "Y" ,"C", "q",  "N", "mc", "mcw", "u",
    "Ht", "av_tax_rate", "T", "I", "B", "BY","TY", "mcww", "G", "τlev", "τprog", "GiniC",
    "GiniX", "sdlogy", "I90share", "I90sharenet", "w90share", "Ygrowth", "Bgrowth",
    "Igrowth", "wgrowth", "Cgrowth", "Tgrowth", "LP", "LPXA", "totRetainedY",
    "union_firm_profits", "unionprofits", "firm_profits", "profits"
]

# ascii names used for cases where unicode doesn't work, e.g., file saves
state_names_ascii = [
    "A", "Z", "ZI", "RB" , "mu","muw", "sigma", "union_Retained", "Retained", "Ylag",
    "Blag", "Tlag", "Ilag", "wlag", "qlag", "Clag", "av_tax_ratelag", "tauproglag",
    "Gshock", "Tprogshock", "Rshock", "Sshock"
]

control_names_ascii = [ 
    "r", "w", "K", "pi" , "piw", "Y" , "C", "q",  "N", "mc", "mcw", "u", "Ht",
    "av_tax_rate", "T", "I", "B", "BY", "TY", "mcww", "G", "taulev", "tauprog",
    "GiniC", "GiniX", "sdlogy", "I90share", "I90sharenet", "w90share", "Ygrowth",
    "Bgrowth", "Igrowth", "wgrowth", "Cgrowth", "Tgrowth", "LP", "LPXA",
    "totRetainedY", "union_firm_profits", "unionprofits", "firm_profits", "profits"
]
# All names in one array
aggr_names          = [state_names; control_names]
aggr_names_ascii    = [state_names_ascii; control_names_ascii]
# Identify distributional summary variables
distr_names         = ["GiniC", "GiniX", "sdlogy", "I90share","I90sharenet",  "w90share"]