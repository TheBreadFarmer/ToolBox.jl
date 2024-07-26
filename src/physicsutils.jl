### Particle Smearing
"""
    smear(particle::LorentzVector, resolution::Float64)

Return smeared reconstructed particle 4-momentum. Resolution is the % resolution in decimal form.
"""
function smear(particle::LorentzVector, resolution::Float64=0.005)
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