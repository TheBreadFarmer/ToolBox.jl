"""
    W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn=0.93891875433, Eb=0.0315)

Compute the hadronic invariant mass W using the approximation that the nucleon remains at rest.

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
function W_approx(l_in::LorentzVector, l_out::LorentzVector; Mn::Float64=AVG_NUCLEON_MASS, Eb::Float64=BINDING_ENERGY)
    # calculate W using approximation. if the argument of sqrt() is negative, return 0.0 instead.
    # `mass(l_in - l_out)^2` is equivalent to Q^2.
    W = try sqrt(-mass(l_in - l_out)^2 + 2*(Mn-Eb)*(energy(l_in) - energy(l_out)) + (Mn-Eb)^2) catch; 0.0 end # units of GeV

    return W
end
# TODO find Eb calls in code and replace with default value.

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
Compute cosine of the angle between two vectors.
"""
function cosbetween(v::LorentzVector,u::LorentzVector)
    v′ = SpatialVector(v)
    u′ = SpatialVector(u)
    return dot(v′, u′) / (norm(v′) * norm(u′))
end