/**
 * \file shift_scoring.hh
 *
 * \brief ShiftScoring — the shift penalty w_s of the bi-alignment model.
 *
 * Copyright (C) Maria Waldl <code@waldl.org>
 */

#ifndef RNASHIFTALIGN_SHIFT_SCORING_HH
#define RNASHIFTALIGN_SHIFT_SCORING_HH

#include <LocARNA/scoring.hh>

namespace RNAShiftAlign {

/**
 * @brief Shift penalty calculator for bi-alignment
 *
 * Implements the shift penalty function w_s from Equation 10 in the paper.
 * The penalty depends on mismatches between gap patterns in the sequence
 * alignment U and structure alignment V.
 *
 * For columns c_U = (c1_U, c2_U) and c_V = (c1_V, c2_V), where each component
 * is 0 (gap) or 1 (content), the shift penalty is:
 *
 *   w_s = Δ × (|c1_U - c1_V| + |c2_U - c2_V|)
 *
 * This gives:
 * - w_s = 0   if columns match exactly (0 mismatches)
 * - w_s = Δ   if one position differs (1 mismatch)
 * - w_s = 2Δ  if both positions differ (2 mismatches)
 *
 * This captures the cost of "shifting" between sequence and structure
 * alignments when their gap patterns are incongruent.
 */
class ShiftScoring {
public:
    using score_t = LocARNA::infty_score_t;

    /**
     * @brief Construct shift scoring with given penalty parameter
     *
     * @param delta Shift penalty parameter (Δ from the paper)
     */
    explicit ShiftScoring(score_t delta);

    /**
     * @brief Compute shift penalty for column pair
     *
     * Computes w_s = Δ × (|c1_U - c1_V| + |c2_U - c2_V|)
     *
     * @param c1_U First component of U column (0 or 1)
     * @param c2_U Second component of U column (0 or 1)
     * @param c1_V First component of V column (0 or 1)
     * @param c2_V Second component of V column (0 or 1)
     *
     * @return Shift penalty w_s ∈ {0, Δ, 2Δ}
     */
    score_t shift_penalty(int c1_U, int c2_U, int c1_V, int c2_V) const;

    /**
     * @brief Get the shift penalty parameter
     *
     * @return Delta value (Δ)
     */
    score_t get_delta() const { return delta_; }

private:
    score_t delta_;  ///< Shift penalty parameter (Δ)
};

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFT_SCORING_HH
