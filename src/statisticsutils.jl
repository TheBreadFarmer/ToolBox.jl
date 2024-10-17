# TODO add full docstrings to all functions


# TODO: add `nothing` bound functionality.
# TODO: consider renaming to something more clear.
# cuts from the lower edge of x1 and the upper edge of x2. i.e. cut does not include the bin starting at x2 or the bin ending at x1.
"""
    efficiency(sample, signal_index, x1, x2; print_output)

Compute the efficiency of selecting signal events from a larger sample for the given cut.
"""
function efficiency(sample::Vector{Hist1D{T}}, signal_index::Int,
                    x1::Real, x2::Real; print_output::Bool=true) where {T}

    # get hist of selected signal events, where bounds are restricted by cut.
    restricted_signal_hist = restrict(sample[signal_index], x1, x2)

    # calculate the number of signal events selected by the cut.
    nevents_signal_cut = integral(restricted_signal_hist)

    # calculate the total number of signal events.
    nevents_signal_total = integral(sample[signal_index])

    # compute efficiency
    _efficiency = nevents_signal_cut / nevents_signal_total

    if print_output == true
        println("efficiency = $(round(_efficiency; digits=4)*100)%")
    end

    # calculate and return the efficiency of these bounds.
    return _efficiency
end

# TODO: add `nothing` bound functionality.
# TODO: consider renaming to something more clear.
# cuts from the lower edge of x1 and the upper edge of x2. i.e. cut does not include the bin starting at x2 or the bin ending at x1.
"""
    purity(sample, signal_index, x1, x2; print_output)

Compute the purity of selected events from given sample events with respect to given signal events.
"""
function purity(sample::Vector{Hist1D{T}}, signal_index::Int, x1::Real, x2::Real;
                print_output::Bool=true) where {T}

    # get hist of selected signal events, where bounds are restricted by cut.
    restricted_signal_hist = restrict(sample[signal_index], x1, x2)

    # calculate the number of signal events selected by the cut.
    nevents_signal_cut = integral(restricted_signal_hist)

    # calculate the total number of events selected by the cut.
    # HOW?
    # 1.create a vector of the number of entries for each restricted hist from `sample`, restricted to select values between cut.
    # 2.sum all elements of that vector to compute the total number of events selected by the cut.
    nevents_all_cut = sum([integral(restrict(sample[i], x1, x2)) for i in eachindex(sample)])

    # compute the purity for these bounds.
    purity = nevents_signal_cut / nevents_all_cut

    if print_output == true
        println("Purity = $(round(purity; digits=4)*100)%")
    end

    # if purity is NaN replace with zero and return, else return purity.
    isnan(purity) && return 0.0
    return purity
end

# FIXME when one of the histograms in the vector is the product of adding two histograms together, the function will fail. investigate and fix.
# TODO URGENT!! DOCUMENT THIS FUNCTION
"""
    maxpurity(sample, signal_index; return_edges, print_output)

Compute maximum possible purity of the given sample with respect to the given signal.
"""
function maxpurity(sample::Vector{Hist1D{T}}, signal_index::Int;
                    return_edges::Bool=false, print_output::Bool=true) where {T}

    highest_purity_value = 0.0

    highest_purity_bin_edges = (0.0, 0.0)

    _edges = binedges(sample[signal_index])[begin:end]

    # creates pairs up to end of range. is itterated in reference to the lower edge.
    each_bin_edges = tuple.(_edges[1:end-1], _edges[2:end])

    for bin in each_bin_edges
        bin_purity = purity(sample, signal_index, bin[1], bin[2]; print_output=false)

        bin_purity > highest_purity_value && 
        (highest_purity_value=bin_purity; highest_purity_bin_edges=bin)
    end

    if print_output == true
        println("Max Purity = $(round(highest_purity_value; digits=4)*100)%")

        if return_edges == false
            return highest_purity_value
        else
            println("    @ Bounds = $(highest_purity_bin_edges)")

            return highest_purity_value, highest_purity_bin_edges
        end
        
    end

    if return_edges == false
        return highest_purity_value
    else
        return highest_purity_value, highest_purity_bin_edges
    end
end

# Optionally supply a lower and upper bound as fourth and fith argument to preview stats within that bound.
# If no bounds are passed, it will compute the best bounds internaly.
"""
    cut_stats(sample, signal_index, target_purity; print_output)

Returns useful stats regarding cuting a sample by W with respect to a given signal.
"""
function cut_stats(sample::Vector{Hist1D{T}}, signal_index::Int,
                   target_purity::Real; print_output::Bool=true) where {T}

    _max = maxpurity(sample, signal_index; return_edges=false, print_output=false)

    target_purity > _max && error("Keyword argument \"target_purity\" must 
        be below the max purity for the given sample and signal.") # TODO break this out into a seperate function and add to any function where a purity selection is required.

    purity, bounds = best_cut_purity(sample, signal_index, target_purity; print_output=false)

    eff = efficiency(sample, signal_index, bounds...; print_output=false)

    width = bounds[2] - bounds[1]
    
    if print_output == true
        println("Purity = $(round(purity*100; digits=4))%")
        println("Efficiency = $(round(eff*100; digits=4))%")
        println("Bounds = $(bounds)")
        println("Width = $(width)")
    end

    return (purity = purity, efficiency = eff, bounds = bounds, width = width)
end

# TODO change function name to optimizecut()
# This function does not account for overflow bins.
"""
    best_cut_purity(sample, signal_index, target_purity; print_output)
"""
function best_cut_purity(sample::Vector{Hist1D{T}}, signal_index::Int,
                         target_purity::Real; print_output::Bool=true) where {T}

    # define variable for chosen purity.
    chosen_purity = 0

    # define variable for best bounds.
    best_bounds = (0,0)

    # get bin stepsize and upper bin edge.
    step_size = step(binedges(sample[signal_index]).uniform_edges)
    upper_edge = last(binedges(sample[signal_index]))-step_size

    # get range for possible lower bound valeus from first bin to last bin.
    lower_bounds = (first(binedges(sample[signal_index])) + step_size):(step_size):(upper_edge)

    # step through each lower bound.
    for lower_bound in lower_bounds # for each lower bound...
        # get range for possible upper bound values from lower_bound to upper_edge.
        upper_bounds = (lower_bound + step_size):(step_size):(upper_edge)

        # step through each upper bound.
        for upper_bound in upper_bounds # for each lower and upper bound pair...
            # compute purity for given upper_bound and lower_bound.
            bound_purity = purity(sample, signal_index, lower_bound, upper_bound; print_output=false)

            # check purity is greater than minimum set by target_purity.
            if bound_purity >= target_purity # purity greater than or equal to target_purity.
                # check width between lower_bound and upper_bound is greater than width between "best" bounds.
                if upper_bound - lower_bound > best_bounds[2] - best_bounds[1]
                    # set best_bounds to lower and upper bounds.
                    best_bounds = (lower_bound, upper_bound)

                    # set chosen_purity to purity.
                    chosen_purity = bound_purity
                end
            end
        end
    end

    # compute width of bounds
    width = best_bounds[2] - best_bounds[1]

    if print_output == true
        println("Purity = $(round(chosen_purity; digits=4))")
        println("Bounds = $(best_bounds)")
        println("Width = $(width)")
    end

    return chosen_purity, best_bounds
end