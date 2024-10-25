
# FIXME check if BINDING_ENERGY is specifically for 12-Carbon.
"Average binding energy for nucleon in (12-Carbon?) nucleus. Units [GeV]"
BINDING_ENERGY = -0.0315 # GeV

"Units [GeV]"
NEUTRON_MASS = 0.9395654205 # GeV
"Units [GeV]"
PROTON_MASS = 0.93827208816 # GeV

"truncated to 10 significant digits to keep accuracy of masses used to compute value. Units [GeV]"
AVG_NUCLEON_MASS = trunc((PROTON_MASS+NEUTRON_MASS)/2,sigdigits=10) # GeV