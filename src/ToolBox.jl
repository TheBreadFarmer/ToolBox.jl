module ToolBox

using LorentzVectors, LorentzVectorHEP, GLMakie, FHist, Distributions


include("./statutils.jl")
include("./plottingutils.jl")


export W_approx, W_exact
export efficiency, purity, max_purity, cut_stats
export smear, get_acceptance
export make_stackedhist





"""
    W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn=0.93891875433, Eb=0.0315)

Compute the approximate hadronic invariant mass W. Return in units passed

```math
W = √{ -(-q)² + 2⋅(Mn-Eb)⋅ω + (Mn-Eb)² }
 q = l_in - l_out
 ω = El_in - El_out
```

W is calculated from the incoming lepton `l_in` and the outgoing 
final state primary lepton `l_out`. 
It is also assumed that the struck nucleon is at rest when struck. 
`Mn` is the mass (in GeV) of the struck nucleon as the average between the 
masses of the neutron and proton reported by PDG. `Eb` is the binding energy term which 
is the energy (in GeV) required to kick a nucleon out of the nucleus (defaults to 0.0).

See also [`W_exact`](@ref)
"""
function W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn::Float64=0.93891875433, Eb::Float64=0.0)
    # calculate W using approximation. if the argument of sqrt() is negative, return 0.0 instead.
    # `mass(l_in - l_out)^2` is equivalent to Q^2.
    W = try sqrt(-mass(l_in - l_out)^2 + 2*(Mn-Eb)*(energy(l_in) - energy(l_out)) + (Mn-Eb)^2) catch; 0.0 end # units of GeV

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
