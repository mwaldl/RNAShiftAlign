#ifndef RNASHIFTALIGN_SHIFT_ALIGNER_HH
#define RNASHIFTALIGN_SHIFT_ALIGNER_HH

#include <string>
#include <memory>

#include <LocARNA/rna_data.hh>
#include <LocARNA/scoring.hh>
#include <LocARNA/aligner_params.hh>
#include <LocARNA/arc_matches.hh>

#include "shiftmatrix_m.hh"
#include "shiftmatrix_d.hh"
#include "shift_scoring.hh"

namespace RNAShiftAlign {

/**
 * @brief Parameters for RNA structure preprocessing
 *
 * Controls base pair probability computation and arc match filtering.
 */
struct PreprocessingParams {
    double min_prob;              ///< Minimum base pair probability (default: 0.0005)
    double max_bps_length_ratio;  ///< Max base pairs / sequence length (default: 0.0)
    int max_diff_am;              ///< Max arc length difference (default: -1 = off)
    int max_diff_at_am;           ///< Max position difference at arc ends (default: -1 = off)

    PreprocessingParams()
        : min_prob(0.0005),
          max_bps_length_ratio(0.0),
          max_diff_am(-1),
          max_diff_at_am(-1) {}
};

/**
 * @brief Parameters for alignment scoring and algorithm
 *
 * Defines match/mismatch scores, gap costs, structure weights, and shift parameters.
 */
struct ScoringParams {
    // Base scoring
    int match;           ///< Base match score (default: 50)
    int mismatch;        ///< Base mismatch score (default: 0)
    int indel;           ///< Gap penalty (default: -150)
    int indel_opening;   ///< Gap opening penalty (default: -500, unused for linear gaps)

    // Structure scoring
    int struct_weight;   ///< Structure contribution weight (default: 200)
    int tau_factor;      ///< Sequence contribution to arc matches (default: 0)

    // Shift parameters
    int delta;           ///< Shift penalty (Δ) for gap pattern mismatches (default: 100)
    int max_shifts;      ///< Maximum allowed shifts (δ_max heuristic) (default: 5)

    ScoringParams()
        : match(50),
          mismatch(0),
          indel(-100),
          indel_opening(-500),
          struct_weight(200),
          tau_factor(0),
          delta(-150),
          max_shifts(5) {}
};

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
     * Performs complete preprocessing pipeline:
     * 1. Computes base pairing probabilities via ViennaRNA partition function
     * 2. Extracts significant base pairs (filtered by min_prob)
     * 3. Computes valid arc matches between sequences
     * 4. Initializes scoring scheme for base matches and arc matches
     * 5. Allocates shift-aware dynamic programming matrices
     *
     * @param seqA First RNA sequence (ACGU alphabet)
     * @param seqB Second RNA sequence (ACGU alphabet)
     * @param prep_params Preprocessing parameters (base pair probabilities, filtering)
     * @param scoring_params Scoring parameters (match/mismatch, gaps, structure, shifts)
     */
    ShiftAligner(const std::string &seqA,
                 const std::string &seqB,
                 const PreprocessingParams &prep_params = PreprocessingParams(),
                 const ScoringParams &scoring_params = ScoringParams());

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
     * @brief Get RNA data for sequence A (base pair probabilities)
     */
    const LocARNA::RnaData& get_rna_dataA() const { return *rna_dataA_; }

    /**
     * @brief Get RNA data for sequence B (base pair probabilities)
     */
    const LocARNA::RnaData& get_rna_dataB() const { return *rna_dataB_; }

    /**
     * @brief Get arc matches between sequences
     */
    const LocARNA::ArcMatches& get_arc_matches() const { return *arc_matches_; }

    /**
     * @brief Get LocARNA scoring object
     */
    const LocARNA::Scoring& get_scoring() const { return *locarna_scoring_; }

    /**
     * @brief Get computed alignment score
     *
     * Only valid after align() has been called.
     */
    score_t get_score() const { return alignment_score_; }

    /**
     * @brief Get alignment U (sequence layer)
     *
     * Only valid after traceback() has been called.
     * @return Pair of aligned sequences (seqA with gaps, seqB with gaps)
     */
    std::pair<std::string, std::string> get_alignment_U() const {
        return {alignment_U_seqA_, alignment_U_seqB_};
    }

    /**
     * @brief Get alignment V (structure layer)
     *
     * Only valid after traceback() has been called.
     * @return Pair of aligned sequences (seqA with gaps, seqB with gaps)
     */
    std::pair<std::string, std::string> get_alignment_V() const {
        return {alignment_V_seqA_, alignment_V_seqB_};
    }

    /**
     * @brief Get formatted alignment string for visualization(y1,y2,ar,br)
     *
     * Returns human-readable alignment showing:
     * - Sequence alignment (U) with sequence and structure for both RNAs
     * - Structure alignment (V) with sequence and structure for both RNAs
     * - Shift annotation showing where U ≠ V
     *
     * Only valid after traceback() has been called.
     * @return Formatted multi-line string
     */
    std::string format_alignment() const;

private:
    // Sequences
    std::unique_ptr<LocARNA::Sequence> seqA_;
    std::unique_ptr<LocARNA::Sequence> seqB_;

    // Preprocessing data
    std::unique_ptr<LocARNA::RnaData> rna_dataA_;
    std::unique_ptr<LocARNA::RnaData> rna_dataB_;
    std::unique_ptr<LocARNA::ArcMatches> arc_matches_;

    // Scoring
    std::unique_ptr<LocARNA::Scoring> locarna_scoring_;
    ShiftScoring shift_scoring_;

    // Heuristic parameters
    size_type max_shifts_;  ///< Maximum allowed number of shifts (δ_max)

    // Dynamic programming matrices
    std::unique_ptr<LocARNA::ShiftMatrixM<score_t>> M_;  ///< Main DP matrix (4D)
    std::unique_ptr<LocARNA::ShiftMatrixD<score_t>> D_;  ///< Structure DP matrix (6D)

    // Results
    score_t alignment_score_;  ///< Computed optimal score

    // Alignment storage (populated by traceback)
    std::string alignment_U_seqA_;  ///< Sequence A in alignment U (with gaps)
    std::string alignment_U_seqB_;  ///< Sequence B in alignment U (with gaps)
    std::string alignment_V_seqA_;  ///< Sequence A in alignment V (with gaps)
    std::string alignment_V_seqB_;  ///< Sequence B in alignment V (with gaps)

    // ========== Core alignment functions (following LocARNA pattern) ==========

    /**
     * @brief Initialize D matrix
     *
     * Creates D matrix with proper dimensions based on arc matches
     * and max_shifts parameter. Must be called before align_D().
     */
    void init_D();

    /**
     * @brief Fill D matrix for all arc matches
     *
     * Iterates over left arc endpoints (al, bl) in descending order.
     * For each (al, bl) and each valid sequence offset (x1, x2):
     * 1. Calls align_in_arcmatch to fill local M for the arc region
     * 2. Calls fill_D_entries to extract D values
     *
     * Inner arcs are processed before outer arcs.
     */
    void align_D();

    /**
     * @brief Fill M matrix for region inside an arc match
     *
     * Initializes M matrix for starting column (x1, x2, al, bl) and fills
     * all entries for positions inside the arc region.
     *
     * @param al Left endpoint of arc in sequence A (structure layer)
     * @param ar Right endpoint of arc in sequence A (structure layer)
     * @param bl Left endpoint of arc in sequence B (structure layer)
     * @param br Right endpoint of arc in sequence B (structure layer)
     * @param x1 Sequence position at left end in A (can shift from al)
     * @param x2 Sequence position at left end in B (can shift from bl)
     * @param y1 Sequence position at right end in A (can shift from ar)
     * @param y2 Sequence position at right end in B (can shift from br)
     */
    void align_in_arcmatch(size_type al, size_type ar,
                           size_type bl, size_type br,
                           size_type x1, size_type x2,
                           size_type y1, size_type y2);

    /**
     * @brief Core recursion: compute M value at position
     *
     * Optimizes over:
     * - Case 1: All 15 unpaired column types
     * - Case 2: All arc matches with right ends at (y3, y4)
     *
     * @param al Left boundary of current arc region (structure A)
     * @param bl Left boundary of current arc region (structure B)
     * @param x1 Starting sequence position in A (at arc left)
     * @param x2 Starting sequence position in B (at arc left)
     * @param y1 Current sequence position in A (sequence layer)
     * @param y2 Current sequence position in B (sequence layer)
     * @param y3 Current structure position in A (structure layer)
     * @param y4 Current structure position in B (structure layer)
     * @return Optimal score for M(y1, y2, y3, y4)
     */
    score_t align_noex(size_type al, size_type bl,
                       size_type x1, size_type x2,
                       size_type y1, size_type y2,
                       size_type y3, size_type y4);

    /**
     * @brief Extract D matrix entries for arc matches with left ends (al, bl)
     *
     * For each arc match with left ends at (al, bl) and each valid
     * sequence offset combination:
     * D(am, z1, z2, y1, y2) = M(y1, y2, ar-1, br-1) + arcmatch_score(am)
     *
     * @param al Left endpoint in sequence A (structure layer)
     * @param bl Left endpoint in sequence B (structure layer)
     * @param x1 Sequence position at left end in A
     * @param x2 Sequence position at left end in B
     */
    void fill_D_entries(size_type al, size_type bl,
                        size_type x1, size_type x2);

    /**
     * @brief Initialize M matrix boundaries for arc region
     *
     * Sets up M matrix entries at boundaries for region defined by
     * structure positions (al, ar, bl, br) and sequence positions
     * (x1, x2) to (y1, y2).
     *
     * @param al Left endpoint of arc in sequence A (structure layer)
     * @param ar Right endpoint of arc in sequence A (structure layer)
     * @param bl Left endpoint of arc in sequence B (structure layer)
     * @param br Right endpoint of arc in sequence B (structure layer)
     * @param x1 Sequence position at left end in A
     * @param x2 Sequence position at left end in B
     * @param y1 Sequence position at right end in A
     * @param y2 Sequence position at right end in B
     */
    void init_M(size_type al, size_type ar,
                size_type bl, size_type br,
                size_type x1, size_type x2,
                size_type y1, size_type y2);

    /**
     * @brief Check if shift between positions is valid
     *
     * @param seq_pos Sequence layer position
     * @param struct_pos Structure layer position
     * @return true if |seq_pos - struct_pos| <= max_shifts_
     */
    bool valid_shift(size_type seq_pos, size_type struct_pos) const {
        return static_cast<int>(std::abs(static_cast<int>(seq_pos) -
                                         static_cast<int>(struct_pos)))
               <= static_cast<int>(max_shifts_);
    }

    // ========== Helper functions ==========

    /**
     * @brief Compute score for unpaired columns (Case 1)
     *
     * Iterates over all 15 valid column types and returns best score.
     */
    score_t compute_unpaired_score(size_type al, size_type bl,
                                   size_type x1, size_type x2,
                                   size_type y1, size_type y2,
                                   size_type y3, size_type y4);

    /**
     * @brief Compute score for paired columns (Case 2)
     *
     * Iterates over all arc matches with right ends at (y3, y4).
     */
    score_t compute_paired_score(size_type al, size_type bl,
                                 size_type x1, size_type x2,
                                 size_type y1, size_type y2,
                                 size_type y3, size_type y4);

    // ========== Legacy functions (used by current align()) ==========

    /**
     * @brief Fill D matrix for all arc matches (legacy)
     *
     * Iterates over all arc matches and fills D matrix entries.
     * Will be replaced by align_D() when restructuring is complete.
     */
    void fill_D();

    /**
     * @brief Fill M matrix for local region (legacy)
     *
     * Fills M matrix for region defined by arc boundaries.
     * Will be replaced by align_in_arcmatch() when restructuring is complete.
     *
     * @param left_A Left boundary in sequence A
     * @param right_A Right boundary in sequence A
     * @param left_B Left boundary in sequence B
     * @param right_B Right boundary in sequence B
     */
    void fill_M_local(size_type left_A, size_type right_A,
                      size_type left_B, size_type right_B);

    /**
     * @brief Fill M matrix for entire sequences with structure (legacy)
     *
     * Top-level alignment call for the full sequence range.
     * Will be replaced by align() calling align_in_arcmatch when restructuring is complete.
     */
    void fill_M_with_structure();
};

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFT_ALIGNER_HH
