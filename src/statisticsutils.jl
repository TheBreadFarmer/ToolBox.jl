"""
    efficiency(sample, signal_index, x1, x2)

Compute the efficiency of selecting signal events from a larger sample for the given cut.

efficiency = selected signal events / total signal events

cuts from the lower edge of x1 and the upper edge of x2. i.e. cut does not include the bin starting at x2 or the bin ending at x1.
"""
function efficiency(sample::Vector{Hist1D{T}}, signal_index::Int, x1::Real, x2::Real) where {T}
    # select signal events within bounds.
    selected_signal_events = restrict(sample[signal_index], x1, x2)

    # compute number of selected signal events.
    nsignal_events_selected = integral(selected_signal_events)

    # compute total number of signal events.
    nsignal_events_total = integral(sample[signal_index])

    # compute efficiency
    _efficiency = nsignal_events_selected / nsignal_events_total

    return _efficiency
end



# TODO add further explanation for what this function does.
"""
    purity(sample, signal_index, x1, x2)

Return purity of selection with respect to signal.

purity = selected signal events / total selected events

TODO add further explanation for what this function does.

cuts from the lower edge of x1 and the upper edge of x2. i.e. cut does not include the bin starting at x2 or the bin ending at x1.
"""
function purity(sample::Vector{Hist1D{T}}, signal_index::Int, x1::Real, x2::Real) where {T}
    # restrict each histogram in sample by bounds x1 and x2.
    selected_events = [restrict(hist, x1, x2) for hist in sample]

    # compute number of total events selected.
    ntotal_events_selected = sum(integral.(selected_events))

    # compute number of signal events selected.
    nsignal_events_selected = integral(selected_events[signal_index])

    # purity = ratio of selected signal to all selected.
    purity = nsignal_events_selected / ntotal_events_selected

    # TODO investigate this issue and fix so that this can be removed and `return purity` can be used instead.
    # if purity == NaN, return 0.0 instead.
    return isnan(purity) ? 0.0 : purity
end



# FIXME when one of the histograms in the vector is the product of adding two histograms together, the function will fail. investigate and fix.
"""
    maxpurity(sample, signal_index) -> max_purity_value, max_purity_bin_edges

Return purity value and edges of bin with highest purity value.
"""
function maxpurity(sample::Vector{Hist1D{T}}, signal_index::Int) where {T}
    max_purity_value = 0.0

    max_purity_bin_edges = (0.0, 0.0)

    sample_bin_edges = binedges(sample[signal_index]).uniform_edges

    lower_edges = sample_bin_edges[1:end-1]
    upper_edges = sample_bin_edges[2:end]

    for (bin_lower_edge,bin_upper_edge) in zip(lower_edges,upper_edges)
        bin_purity = purity(sample, signal_index, bin_lower_edge, bin_upper_edge)

        if bin_purity > max_purity_value
            max_purity_value = bin_purity

            max_purity_bin_edges = (bin_lower_edge, bin_upper_edge)
        end
    end

    return max_purity_value, max_purity_bin_edges
end





# This function does not account for overflow bins.
# TODO make docstring more clear.
"""
    optimal_bounds(sample, signal_index, target_purity)

Return bounds with purity closest (from above) to `target_purity`.
"""
function optimal_bounds(sample::Vector{Hist1D{T}}, signal_index::Int,
                        target_purity::Real) where {T}
    best_bound = (0,0)
    best_bound_purity = 0

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
            bound_purity = purity(sample, signal_index, lower_bound, upper_bound)

            # check purity is greater than minimum set by target_purity.
            if bound_purity >= target_purity # purity greater than or equal to target_purity.
                # check width between lower_bound and upper_bound is greater than width between "best" bounds.
                if upper_bound - lower_bound > best_bound[2] - best_bound[1]
                    # set best_bounds to lower and upper bounds.
                    best_bound = (lower_bound, upper_bound)

                    # set best_bound_purity to purity.
                    best_bound_purity = bound_purity
                end
            end
        end
    end

    return chosen_purity, best_bounds
end


# TODO refactor this function into a printing function.
# Optionally supply a lower and upper bound as fourth and fith argument to preview stats within that bound.
# If no bounds are passed, it will compute the best bounds internaly.
"""
    cut_stats(sample, signal_index, target_purity; print_output)

Returns useful stats regarding cuting a sample by W with respect to a given signal.
"""
function cut_stats(sample::Vector{Hist1D{T}}, signal_index::Int,
                   target_purity::Real; print_output::Bool=true) where {T}

    _max = maxpurity(sample, signal_index)

    target_purity > _max && error("Keyword argument \"target_purity\" must 
        be below the max purity for the given sample and signal.") # TODO break this out into a seperate function and add to any function where a purity selection is required.

    purity, bounds = optimal_bounds(sample, signal_index, target_purity)

    eff = efficiency(sample, signal_index, bounds...)

    width = bounds[2] - bounds[1]
    
    if print_output == true
        println("Purity = $(round(purity*100; digits=4))%")
        println("Efficiency = $(round(eff*100; digits=4))%")
        println("Bounds = $(bounds)")
        println("Width = $(width)")
    end

    return (purity = purity, efficiency = eff, bounds = bounds, width = width)
end