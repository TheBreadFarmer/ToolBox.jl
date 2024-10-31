"""
    W_approx(l_in, l_out; Mn=AVG_NUCLEON_MASS, Eb=BINDING_ENERGY)

Compute approximation of the hadronic invariant mass `W` assuming the struck nucleon is at rest with a correction for binding energy `Eb`.

```math
W = √{ -(-q)² + 2⋅(Mn+Eb)⋅ω + (Mn+Eb)² }
 q = (l_in - l_out)
 ω = El_in - El_out
```

# Arguments
- `l_in::LorentzVector`: the incoming lepton's 4-momentum.
- `l_out::LorentzVector`: the primary final-state outgoing lepton's 4-momentum.
- `Mn::Real=AVG_NUCLEON_MASS`: the mass of the struck nucleon.
- `Eb::Real=BINDING_ENERGY`: the binding energy of the nucleus.


See also [`W_exact`](@ref)
"""
function W_approx(l_in::Vec4, l_out::Vec4; Mn::Real=AVG_NUCLEON_MASS, Eb::Real=BINDING_ENERGY) #FIXME bad function name
    # energy transfered
    ω = energy(l_in) - energy(l_out)
    # transfered four-momentum. `mass(l_in - l_out)` is equivalent to Q.
    Q = mass(l_in - l_out)
    # struck-nucleon mass term with correction for binding energy.
    M = Mn+Eb

    # hadronic invariant mass squared.
    W2 = M^2 - Q^2 + 2*M*ω

    # return 0.0 if W2 is negative.
    return W2 < 0.0 ? 0.0 : sqrt(W2)
end


"""
    W_exact(l_in::LorentzVector, l_out::LorentzVector, n_in::LorentzVector)

Compute the hadronic invariant mass `W` without approximation.

```math
W = √{ -(-q)² + 2( q ⋅ n_in ) + Mn² }
q = l_in - l_out
```

# Arguments
- `l_in::LorentzVector`: the incoming lepton's 4-momentum.
- `l_out::LorentzVector`: the primary final-state outgoing lepton's 4-momentum.
- `n_in::LorentzVector`: the struck nucleon.

See also [`W_approx`](@ref)
"""
function W_exact(l_in::LorentzVector, l_out::LorentzVector, n_in::LorentzVector) #TODO make this more readable.
    # calculate W using exact values. `mass(l_in - l_out)^2` is equivalent to Q^2.
    W = sqrt(-mass(l_in - l_out)^2 + 2*dot((l_in - l_out), n_in) + mass(n_in)^2)


    return W
end




"""
    cosbetween(v, u)

Compute cosine of the angle between two vectors.
"""
function cosbetween(v::SpatialVector,u::SpatialVector)
    return dot(v, u) / (norm(v) * norm(u))
end
cosbetween(v::LorentzVector,u::LorentzVector) = cosbetween(SpatialVector(v),SpatialVector(u))




"""
    predict_pion_energy(lepton_in, lepton_out_smeared, fs_lead_proton, Mn=NEUTRON_MASS, Eb=BINDING_ENERGY)

Compute predicted pion energy acording to Eπ = (Ee - Ee′) + Mn + Eb - Ep, for this interaction: e + n -> e′ + p + π⁻

The incoming nucleon must be a neutron for the interaction to be valid.
"""
function predict_pion_energy(lepton_in, lepton_out_smeared, fs_lead_proton; Mn=NEUTRON_MASS, Eb=BINDING_ENERGY)
    Ee = energy(lepton_in)
    Ee′ = energy(lepton_out_smeared)
    Ep = energy(fs_lead_proton)

    return (Ee - Ee′) + Mn + Eb - Ep
end