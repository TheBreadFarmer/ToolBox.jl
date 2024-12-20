module ToolBox

import Distributions: Normal

using LorentzVectors, LorentzVectorHEP, FHist



export AVG_NUCLEON_MASS, BINDING_ENERGY, NEUTRON_MASS, PROTON_MASS
export W_CUT_BOUNDS, RESIDUAL_CPI_ENERGY_BOUNDS

include("./constants/generalconstants.jl")
include("./constants/physicalconstants.jl")

export efficiency, purity, maxpurity
export optimal_bounds, cut_stats

include("./statisticsutils.jl")

export W_approx, W_exact
export cosbetween, predict_pion_energy

include("./mathutils.jl")

export smear, get_acceptance

include("./physicsutils.jl")

export w_cut, residual_cpi_energy_cut, equalsany, particles, get_lead_particle, reconstruct_particles

include("./misc.jl")

end