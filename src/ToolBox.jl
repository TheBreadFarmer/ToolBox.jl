module ToolBox

using LorentzVectors, LorentzVectorHEP, GLMakie, FHist, Distributions


export W_approx, W_exact
export efficiency, purity, maxpurity, cut_stats
export smear, get_acceptance


include("physicalconstants.jl")
include("statisticsutils.jl")
include("mathutils.jl")
include("physicsutils.jl")

end