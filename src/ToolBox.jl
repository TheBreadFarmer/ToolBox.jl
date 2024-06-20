module ToolBox

using LorentzVectors, LorentzVectorHEP, GLMakie, FHist, Distributions


include("./statutils.jl")
include("./plottingutils.jl")


export W_approx, W_exact
export efficiency, purity, max_purity, cut_stats
export smear, get_acceptance
export make_stackedhist





"""
    W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn::Float64=938.91875433)

Compute the approximate hadronic invariant mass W.

```math
W = √{ -(-q)² + 2⋅Mn⋅ω + Mn² }
 q = l_in - l_out
 ω = El_in - El_out
```

W is calculated from the incoming lepton `l_in` and the outgoing final state primary lepton `l_out`. 
It is also assumed that the struck nucleon is at rest when struck. 
`Mn` is the mass (in MeV) of the struck nucleon as the average between the masses of the neutron and proton reported by PDG.
The units of the returned value for W will be the same as the units of the argument value.

See also [`W_exact`](@ref)
"""
function W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn::Float64=938.91875433)
    # calculate W using approximation. if the argument of sqrt() is negative, return 0.0 instead.
    # `mass(l_in - l_out)^2` is equivalent to Q^2.
    W = try sqrt(-mass(l_in - l_out)^2 + 2*(Mn/1000)*(energy(l_in) - energy(l_out)) + (Mn/1000)^2) catch; 0.0 end # units of GeV


    return W
end

"""
    W_exact(l_in::LorentzVector, l_out::LorentzVector, n_in::LorentzVector)

Compute the hadronic invariant mass W.

```math
W = √{ -(-q)² + 2( q ⋅ n_in ) + Mn² }
q = l_in - l_out
```

W is calculated from the four-momentum of the incoming lepton `l_in`,
the outgoing final state primary lepton `l_out`,
and the struck/incoming nucleon `n_in`.
This function does not assume that the struck nucleon is at rest.
The units of the returned value for W will be the same as the units of the argument value.

See also [`W_approx`](@ref)
"""
function W_exact(l_in::LorentzVector, l_out::LorentzVector, n_in::LorentzVector)
    # calculate W using exact values. `mass(l_in - l_out)^2` is equivalent to Q^2.
    W = sqrt(-mass(l_in - l_out)^2 + 2*dot((l_in - l_out), n_in) + mass(n_in)^2) # units of GeV


    return W
end



"""
    research_dir()

Return absolute path to primary directory of research project. Throw error if no such path exists.

CAUTION: Make sure this path is kept up to date if the directory is moved, altered, or folder names are changed.

See also [`plot_tempdir`](@ref)
"""
function research_dir()
    path = raw"C:/Users/phfit/Desktop/Summer 2024/Research/"
    isdir(path) ? string(path) : error("$(path) could not be found.\nPlease update variable `path` in function $(string(nameof(research_dir)))()\n    @ \"$(functionloc(research_dir)[1])\":$(functionloc(research_dir)[2])")
end


"""
    plot_tempdir()

Return absolute path to "~/plotlib/temp" directory. Throw error if no such path exists.

CAUTION: Make sure this path is kept up to date if the directory is moved, altered, or folder names are changed.

See also [`research_dir`](@ref)
"""
function plot_tempdir()
    path = string(research_dir(), raw"/plotlib/temp/")
    isdir(path) ? string(path) : error("$(path) could not be found.\nPlease update variable `path` in function $(string(nameof(plot_tempdir)))()\n    @ \"$(functionloc(plot_tempdir)[1])\":$(functionloc(plot_tempdir)[2])")
end



### Particle Smearing
"""
    smear(particle::LorentzVector, resolution::Float64)

Return smeared reconstructed particle 4-momentum. Resolution is the % resolution in decimal form.
"""
function smear(particle::LorentzVector, resolution::Float64=0.0005)
    # compute particles true total momentum.
    p_true = LorentzVectorHEP.mag(particle) # units of GeV

    # compute standard deviation from resolution and true total momentum.
    σ = p_true * resolution # units of GeV

    # randomly smear the total momentum.
    p_reco = rand(Normal(p_true, σ)) # units of GeV

    # compute smear ratio.
    r = p_reco / p_true

    # reconstruct smeared particle energy.
    m2 = LorentzVectorHEP.mass2(particle)
    E_reco = sqrt(m2 + r^2*p_true^2)

    # reconstruct smeared particle 3-momentum components.
    x_reco, y_reco, z_reco = r*particle.x, r*particle.y, r*particle.z

    # pack reconstructed components into 4-momentum.
    particle_reco = LorentzVector(E_reco, x_reco, y_reco, z_reco)

    # return reconstructed particle
    return particle_reco
end

"""
    smear(p_true::Float64, resolution::Float64)

Return smeared reconstructed total momentum magnitude. Resolution is the % resolution in decimal form.
"""
function smear(p_true::Float64, resolution::Float64)
    # compute standard deviation from resolution and true total momentum.
    σ = p_true * resolution # units of GeV

    # randomly smear the total momentum.
    p_reco = rand(Normal(p_true, σ)) # units of GeV

    # return the smeared reconstruction of total momentum.
    return p_reco
end



### particle acceptance
function get_acceptance(lv::LorentzVector)
    cos_theta = LorentzVectorHEP.CosTheta(lv)

    return ifelse(0.35 < cos_theta && cos_theta < 0.98, 1.0, 0.0)
end

end
