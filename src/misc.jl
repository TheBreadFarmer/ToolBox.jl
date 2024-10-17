# everything found in here should be moved to a more relevant file with a relevant name for better organization.

"""
    w_cut(w; bounds=W_CUT_BOUNDS) -> Bool

Return `true` if hadronic invariant mass `w` is within the bounds (inclusive) specified by `bounds` which default to `W_CUT_BOUNDS`.
"""
function w_cut(w; bounds=W_CUT_BOUNDS)
    return bounds[1] ≤ w ≤ bounds[2]
end