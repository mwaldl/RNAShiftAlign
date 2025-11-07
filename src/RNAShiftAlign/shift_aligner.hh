#ifndef RNASHIFTALIGN_SHIFT_ALIGNER_HH
#define RNASHIFTALIGN_SHIFT_ALIGNER_HH

#include <string>
#include <memory>

#include <LocARNA/rna_data.hh>
#include <LocARNA/scoring.hh>
#include <LocARNA/aligner_params.hh>

#include "shiftmatrix_m.hh"
#include "shift_scoring.hh"

namespace RNAShiftAlign {

/**
 * @brief Sankoff-style bi-aligner for detecting incongruent RNA evolution
 *
 * Implements the shift-aware bi-alignment algorithm from:
 * "Incongruences Between Sequence and Secondary Structure Alignments"
 *
 * The algorithm simultaneously computes:
 * - U: sequence-driven alignment
 * - V: structure-driven alignment
 * - W: shift tracking between U and V
 *
 * Key features:
 * - 4D dynamic programming with ShiftMatrixM (M matrix)
 * - 6D structure matching with ShiftMatrixD (D matrix, future)
 * - Shift penalty scoring for gap pattern mismatches
 * - δ_max heuristic for O(n²) complexity
 */
class ShiftAligner {
public:
    using score_t = LocARNA::infty_score_t;
    using size_type = size_t;

    /**
     * @brief Construct aligner for two RNA sequences
     *
     * For the skeleton: stores sequences and allocates DP matrices.
     * TODO: Add LocARNA preprocessing when implementing recursion:
     * - Compute base pairing probabilities (RnaData)
     * - Initialize scoring scheme (Scoring with ArcMatches)
     *
     * @param seqA First RNA sequence
     * @param seqB Second RNA sequence
     * @param delta Shift penalty parameter (Δ)
     * @param max_shifts Maximum allowed number of shift positions (δ_max heuristic)
     */
    ShiftAligner(const std::string &seqA,
                 const std::string &seqB,
                 score_t delta,
                 size_type max_shifts);

    /**
     * @brief Perform bi-alignment computation
     *
     * Fills dynamic programming matrices and computes optimal score.
     * Currently implements unpaired-only alignment (no structure).
     *
     * @return Optimal alignment score
     */
    score_t align();

    /**
     * @brief Perform traceback to reconstruct alignments
     *
     * Recovers U, V, and W alignments from DP matrices.
     * To be implemented after forward recursion is complete.
     */
    void traceback();

    /**
     * @brief Get sequence A
     */
    const LocARNA::Sequence& get_seqA() const { return *seqA_; }

    /**
     * @brief Get sequence B
     */
    const LocARNA::Sequence& get_seqB() const { return *seqB_; }

    /**
     * @brief Get shift penalty parameter
     */
    score_t get_delta() const { return shift_scoring_.get_delta(); }

    /**
     * @brief Get maximum number of allowed shifts
     *
     * Returns δ_max, the constraint on |y1-y3| and |y2-y4|.
     */
    size_type get_max_shifts() const { return max_shifts_; }

    /**
     * @brief Get computed alignment score
     *
     * Only valid after align() has been called.
     */
    score_t get_score() const { return alignment_score_; }

private:
    // Sequences
    std::unique_ptr<LocARNA::Sequence> seqA_;
    std::unique_ptr<LocARNA::Sequence> seqB_;

    // Preprocessing data (to be initialized when implementing recursion)
    std::unique_ptr<LocARNA::RnaData> rna_dataA_;
    std::unique_ptr<LocARNA::RnaData> rna_dataB_;

    // Scoring
    std::unique_ptr<LocARNA::Scoring> locarna_scoring_;
    ShiftScoring shift_scoring_;

    // Heuristic parameters
    size_type max_shifts_;  ///< Maximum allowed number of shifts (δ_max)

    // Dynamic programming matrices
    std::unique_ptr<LocARNA::ShiftMatrixM<score_t>> M_;  ///< Main DP matrix

    // Results
    score_t alignment_score_;  ///< Computed optimal score

    /**
     * @brief Fill M matrix for unpaired-only alignment
     *
     * Implements Case 1 of Equation 17 (no base pair matches).
     * Future: Will be extended with Case 2 for structure.
     */
    void fill_M_unpaired();
};

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFT_ALIGNER_HH
