
# FIXME check if BINDING_ENERGY is specifically for 12-Carbon.
"Average binding energy for nucleon in (12-Carbon?) nucleus. Units [GeV]"
const BINDING_ENERGY = -0.0315 # GeV

"Units [GeV]"
const NEUTRON_MASS = 0.9395654205 # GeV
"Units [GeV]"
const PROTON_MASS = 0.93827208816 # GeV

"truncated to 10 significant digits to keep accuracy of masses used to compute value. Units [GeV]"
const AVG_NUCLEON_MASS = trunc((PROTON_MASS+NEUTRON_MASS)/2,sigdigits=10) # GeV