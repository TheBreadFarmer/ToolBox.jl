module ToolBox

using LorentzVectors, LorentzVectorHEP, GLMakie, FHist, Distributions


export W_approx, W_exact
export efficiency, purity, max_purity, cut_stats
export smear, get_acceptance



include("statisticsutils.jl")
include("mathutils.jl")
include("physicsutils.jl")











end
