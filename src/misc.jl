# everything found in here should be moved to a more relevant file with a relevant name for better organization.

"""
    w_cut(w; bounds=W_CUT_BOUNDS) -> Bool

Return `true` if hadronic invariant mass `w` is within the bounds (inclusive) specified by `bounds`.

Values for the hadronic invariant mass `w` that are within the bounds are considered passing.
"""
function w_cut(w; bounds=W_CUT_BOUNDS)
    return bounds[1] ≤ w ≤ bounds[2]
end


"""
    residual_cpi_energy_cut(ΔE; bounds=RESIDUAL_CPI_ENERGY_BOUNDS) -> Bool

Return `true` if residual charged-pion energy `ΔE` is within the bounds (inclusive) specified by `bounds`.
"""
function residual_cpi_energy_cut(ΔE; bounds=RESIDUAL_CPI_ENERGY_BOUNDS)
    return bounds[1] ≤ ΔE ≤ bounds[2]
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
    get_lead_particle(LVs, weights)

Return particle LV from LVs with highest magnitude.
"""
function get_lead_particle(LVs::Vector{LorentzVector{T}}, weights) where T <: Number
    # compute the momentum of every LV in LVs
    momenta = LorentzVectorHEP.mag.(LVs)

    # get just the index of the max momentum magnitude.
    _, i = findmax(momenta)

    # return LorentzVector from LVs for that index.
    return LVs[i], weights[i]
end


"""
    reconstruct_particles(MCcodes::Vector{Int}, pdg, wgt, E, px, py, pz)

Reconstruct only particles of specified MonteCarlo ID code, given its weight passes reconstruction test.

TODO: Further explain the reconsstruction test process.
"""
function reconstruct_particles(MCcodes::Vector{Int}, pdg, wgt, E, px, py, pz)
    Base.depwarn("Replace call to `reconstruct_particles()` with the two argument 
    method `reconstruct_particle(MCcodes, event)`!", :reconstruct_particles; force=true)
    # create empty vector for reconstructed particles.
    recoed_particles = LorentzVector{Float64}[]

    # get index in pdg code list of each particle equal to code.
    particle_indices = findall(x->(equalsany(MCcodes, x)), pdg)

    # check weight against random number for each element in particle_indices.
    for i in particle_indices
        wgt[i] < rand() && continue # if rand is greater, skip.
        # if wgt was greater, reconstruct particle.
        push!(recoed_particles, LorentzVector(E[i], px[i],py[i],pz[i]))
    end

    return recoed_particles
end

"""
    reconstruct_particles(MCcodes::Vector{Int}, event::LazyEvent)

Reconstruct only particles of specified MonteCarlo ID code, given its weight passes reconstruction test.

Also return the weights for reconstructed particles.

TODO: Further explain the reconstruction test process.
"""
function reconstruct_particles(MCcodes::Vector{Int}, event)
    # create empty vector for reconstructed particles.
    recoed_particles = LorentzVector{Float64}[]
    recoed_weights = Float64[]

    # get index in pdgf code list of each particle equal to a code in MCcodes.
    particle_indices = findall(x->(equalsany(MCcodes, x)), event.pdgf)

    # check weight against random number for each element in particle_indices.
    (; wgtf, Ef, pxf, pyf, pzf) = event
    for i in particle_indices
        wgtf[i] < rand() && continue # if rand is greater, skip.
        # if wgtf was greater, reconstruct particle.
        push!(recoed_particles, LorentzVector(Ef[i],pxf[i],pyf[i],pzf[i]))
        push!(recoed_weights, wgtf[i])
    end

    return recoed_particles, recoed_weights
end