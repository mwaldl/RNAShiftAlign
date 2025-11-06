#ifndef RNASHIFTALIGN_SHIFT_SCORING_HH
#define RNASHIFTALIGN_SHIFT_SCORING_HH

#include <LocARNA/scoring.hh>

namespace RNAShiftAlign {

/**
 * @brief Column types for bi-alignment gap patterns
 *
 * Represents the five possible column types in each alignment (U or V)
 * as defined in Equation 10 of the paper:
 * - MATCH: (i,j) or (k,l) - both sequences contribute
 * - DEL_A: (i,-) or (k,-) - gap in second sequence
 * - INS_A: (-,j) or (-,l) - gap in first sequence
 * - DEL_B: Same as DEL_A (for symmetry in naming)
 * - INS_B: Same as INS_A (for symmetry in naming)
 */
enum class ColumnType {
    MATCH,   ///< Both sequences aligned: (i,j) or (k,l)
    DEL_A,   ///< Gap in sequence A: (i,-) or (k,-)
    INS_A,   ///< Gap in sequence A: (-,j) or (-,l)
    DEL_B,   ///< Gap in sequence B: (i,-) or (k,-) [same as DEL_A]
    INS_B    ///< Gap in sequence B: (-,j) or (-,l) [same as INS_A]
};

/**
 * @brief Shift penalty calculator for bi-alignment
 *
 * Implements the shift penalty function w_s from Equation 10 in the paper.
 * The penalty depends on mismatches between gap patterns in the sequence
 * alignment U and structure alignment V:
 *
 * - w_s = 0   if c_U == c_V (gap patterns match)
 * - w_s = Δ   if one is MATCH and the other is a gap
 * - w_s = 2Δ  if both are gaps but in different sequences
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
     * @brief Compute shift penalty for column type pair
     *
     * Implements the 15×15 lookup table from Equation 10.
     *
     * @param c_U Column type in sequence alignment U
     * @param c_V Column type in structure alignment V
     *
     * @return Shift penalty w_s(c_U, c_V) ∈ {0, Δ, 2Δ}
     */
    score_t shift_penalty(ColumnType c_U, ColumnType c_V) const;

    /**
     * @brief Get the shift penalty parameter
     *
     * @return Delta value (Δ)
     */
    score_t get_delta() const { return delta_; }

private:
    score_t delta_;  ///< Shift penalty parameter (Δ)

    /**
     * @brief Check if column type is a gap (not MATCH)
     *
     * @param c Column type to check
     * @return true if c is DEL_A, INS_A, DEL_B, or INS_B
     */
    static bool is_gap(ColumnType c);

    /**
     * @brief Check if two gap types match
     *
     * Returns true if both gaps are in the same sequence
     * (both DEL or both INS).
     *
     * @param c_U Column type in U (must be a gap)
     * @param c_V Column type in V (must be a gap)
     * @return true if gap types match
     */
    static bool gap_types_match(ColumnType c_U, ColumnType c_V);
};

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFT_SCORING_HH
