
BINDING_ENERGY = -0.0315 # GeV

NEUTRON_MASS = 0.9395654205 # GeV
PROTON_MASS = 0.93827208816 # GeV

# truncated to 10 significant digits to keep accuracy of masses used to compute value.
AVG_NUCLEON_MASS = trunc((PROTON_MASS+NEUTRON_MASS)/2,sigdigits=10) # GeV