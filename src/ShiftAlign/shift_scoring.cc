/**
 * \file shift_scoring.cc
 *
 * \brief Implementation of the shift penalty w_s.
 *
 * Copyright (C) Maria Waldl <code@waldl.org>
 */

#include "shift_scoring.hh"

namespace RNAShiftAlign {

ShiftScoring::ShiftScoring(score_t delta)
    : delta_(delta) {
}

ShiftScoring::score_t
ShiftScoring::shift_penalty(int c1_U, int c2_U, int c1_V, int c2_V) const {
    // Manhattan distance between the two layers' gap patterns, scaled by Δ:
    //   w_s = Δ × (|c1_U - c1_V| + |c2_U - c2_V|)
    //
    // Examples:
    //   (1,1) vs (1,1) → Δ × (0 + 0) = 0   (match)
    //   (1,1) vs (1,0) → Δ × (0 + 1) = Δ   (one mismatch)
    //   (1,1) vs (0,0) → Δ × (1 + 1) = 2Δ  (two mismatches)
    //
    // Compute penalty as delta_ added once per mismatch in A (c1) and mismatch in B (c2)
    score_t penalty(0);
    if (c1_U != c1_V) penalty = penalty + delta_;
    if (c2_U != c2_V) penalty = penalty + delta_;
    return penalty;
}

} // namespace RNAShiftAlign