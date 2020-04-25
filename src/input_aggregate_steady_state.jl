# Set aggregate steady state variabel values
ASS       = 1.0
ZSS       = 1.0
ZISS      = 1.0
μSS       = m_par.μ
μwSS      = m_par.μw
τprogSS       = m_par.τ_prog
τlevSS       = m_par.τ_lev

σSS       = 1.0
τprog_obsSS   = 1.0
GshockSS  = 1.0
RshockSS  = 1.0
TlevshockSS  = 1.0
TprogshockSS  = 1.0

SshockSS  = 1.0
rSS       = 1.0 + interest(KSS,1.0 / m_par.μ, NSS, m_par)
πSS       = 1.0
πwSS      = 1.0
CSS       = (YSS - m_par.δ_0*KSS - GSS)
qSS       = 1.0
mcSS      = 1.0 ./ m_par.μ
mcwSS     = 1.0 ./ m_par.μw
mcwwSS    = wSS*mcwSS
uSS       = 1.0
profitsSS = ProfitsSS

BYSS   = BSS/YSS

TYSS   = TSS/YSS
TlagSS = TSS

YlagSS = YSS
BlagSS = BSS
ZlagSS = ZSS
GlagSS = GSS
IlagSS = ISS
wlagSS = wSS
qlagSS = qSS
NlagSS = NSS
ClagSS = CSS
πlagSS = πSS
σlagSS = σSS
rlagSS = rSS
RBlagSS = RBSS
av_tax_ratelagSS = av_tax_rateSS
τproglagSS = τprogSS

mcwwlagSS= mcwwSS
YgrowthSS = 1.0
BgrowthSS = 1.0
ZgrowthSS = 1.0
GgrowthSS = 1.0
IgrowthSS = 1.0
wgrowthSS = 1.0
qgrowthSS = 1.0
NgrowthSS = 1.0
CgrowthSS = 1.0
πgrowthSS = 1.0
σgrowthSS = 1.0
τproggrowthSS = 1.0
rgrowthSS = 1.0
RBgrowthSS = 1.0
mcwwgrowthSS= 1.0
TgrowthSS = 1.0
HtSS = 1.0
