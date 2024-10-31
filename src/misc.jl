# everything found in here should be moved to a more relevant file with a relevant name for better organization.

"""
    w_cut(w; bounds=W_CUT_BOUNDS) -> Bool

Return `true` if hadronic invariant mass `w` is within the bounds (inclusive) specified by `bounds` which default to `W_CUT_BOUNDS`.
"""
function w_cut(w; bounds=W_CUT_BOUNDS)
    return bounds[1] ≤ w ≤ bounds[2]
end


"""
    equalsany(itr,val)

Returns `true` if any value of `itr` is equal to `val` acording to `==`.
"""
function equalsany(itr, val)
    in(1, val .== itr)
end



"""
    particles(codes, pdg, E, px, py, pz)

Return vector of LVs for particles with `pdg` code equal to any value in `codes`.

# Arguments
- `codes::`
"""
function particles(codes, pdg, E, px, py, pz)
    # create empty vector for particles.
    particles = LorentzVector{Float64}[]

    # get index in lists of each particle with matching code.
    particle_indices = findall(x->(equalsany(codes, x)), pdg)

    # for each particle in particle_indices, create LV and push.
    for i in particle_indices
        push!(particles, LorentzVector(E[i],px[i],py[i],pz[i]))
    end

    return particles
end


"""
    get_lead_particle(LVs)

Return particle LV from LVs with highest magnitude.
"""
function get_lead_particle(LVs::Vector{LorentzVector{T}}) where T <: Number
    # compute the momentum of every LV in LVs
    momenta = LorentzVectorHEP.mag.(LVs)

    # get just the index of the max momentum magnitude.
    _, i = findmax(momenta)

    # return LorentzVector from LVs for that index.
    return LVs[i]
end


"""
    reconstruct_particles(code::Int, pdg, wgt, E, px, py, pz)

Reconstruct only particles with pdg code equal to code. acording to a given weight.
"""
function reconstruct_particles(codes, pdg, wgt, E, px, py, pz)
    # create empty vector for reconstructed particles.
    recoed_particles = LorentzVector{Float64}[]

    # get index in pdg code list of each particle equal to code.
    particle_indices = findall(x->(equalsany(codes, x)), pdg)

    # check weight against random number for each element in particle_indices.
    for i in particle_indices
        wgt[i] < rand() && continue # if rand is greater, skip.
        # if wgt was greater, reconstruct particle.
        push!(recoed_particles, LorentzVector(E[i], px[i],py[i],pz[i]))
    end

    return recoed_particles
end