/**
 * \file shift_aligner.hh
 *
 * \brief ShiftAligner — Sankoff-style RNA bi-aligner, plus its parameter structs.
 *
 * Declares PreprocessingParams and ScoringParams (which hold the tool's default
 * parameter values) and the ShiftAligner class that performs the bi-alignment.
 *
 * Copyright (C) Maria Waldl <maria@bioinf.uni-leipzig.de>
 */

#ifndef RNASHIFTALIGN_SHIFT_ALIGNER_HH
#define RNASHIFTALIGN_SHIFT_ALIGNER_HH

#include <string>
#include <memory>
#include <iostream>
#include <ctime>

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
 *
 * @note These defaults are the single source of truth: they are exactly the
 *       values documented in `rnashiftalign --help`. The CLI only overrides a
 *       field when the corresponding option is given.
 */
struct PreprocessingParams {
    double min_prob;              ///< Minimum base pair probability (default: 0.001)
    double max_bps_length_ratio;  ///< Max base pairs / sequence length (default: 0.0 = off)
    int max_diff_am;              ///< Max arc length difference (default: 100; -1 = off)
    int max_diff_at_am;           ///< Max position difference at arc ends (default: -1 = off)

    PreprocessingParams()
        : min_prob(0.001),
          max_bps_length_ratio(0.0),
          max_diff_am(100),
          max_diff_at_am(-1) {}
};

/**
 * @brief Parameters for alignment scoring and algorithm
 *
 * Defines match/mismatch scores, gap costs, structure weights, and shift parameters.
 *
 * @note Defaults scoring values.
 *       The CLI only overrides a
 *       field when the corresponding option is given.
 */
struct ScoringParams {
    // Base scoring
    int match;           ///< Base match score (default: 50)
    int mismatch;        ///< Base mismatch score (default: 0)
    int indel;           ///< Gap penalty (default: -100)
    int indel_opening;   ///< Gap opening penalty (default: -500, not used in RNAshiftalign)

    // Structure scoring
    int struct_weight;   ///< Structure contribution weight (default: 150)
    int tau_factor;      ///< Sequence contribution to arc matches (0 in RNAsgiftalign)

    // Shift parameters
    int delta;           ///< Shift penalty (Δ) per differing gap-pattern position (default: -200)
    int max_shifts;      ///< Maximum allowed shift (δ_max heuristic) (default: 1)

    ScoringParams()
        : match(50),
          mismatch(0),
          indel(-100),
          indel_opening(-500), // not used in RNAshiftalign, no affine gap costs implemented
          struct_weight(150),
          tau_factor(0),      // 0 in RNAshiftalign, sequence should not contribute to structure score
          delta(-200),
          max_shifts(1) {}
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
 * - W: per-column shift scoring, linking U and V; penalzing columns where U and V disagree
 *
 * A column advances four cursors, one per layer and sequence:
 * (c1, c2) in the sequence layer and (c3, c4) in the structure layer, each 0
 * (gap) or 1 (consume). All 15 combinations except (0,0,0,0) are valid. A column
 * whose two layers disagree is a *shift* and is charged the penalty w_s.
 *
 * Structure and sequence cursors may drift apart by at most δ_max
 * (@ref ScoringParams::max_shifts), which bounds the matrix sizes.
 *
 * Two DP matrices:
 * - ShiftMatrixM (4D): best score for the aligned prefixes inside one arc match region
 * - ShiftMatrixD (6D): best score for the contents of a matched arc pair
 *
 * Typical use:
 * @code
 *   ShiftAligner aligner(seqA, seqB, prep_params, scoring_params);
 *   auto score = aligner.align();     // forward pass
 *   aligner.traceback();              // reconstruct U, V, consensus structure
 *   std::cout << aligner.format_alignment();
 * @endcode
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
     * @param verbose Verbosity level (0=silent, 1=basic, 2=detailed, 3=debug)
     */
    ShiftAligner(const std::string &seqA,
                 const std::string &seqB,
                 const PreprocessingParams &prep_params = PreprocessingParams(),
                 const ScoringParams &scoring_params = ScoringParams(),
                 int verbose = 0);

    /**
     * @brief Run the forward pass and return the optimal score
     *
     * Two phases, following the Sankoff ordering:
     * 1. align_D() fills D for every arc match, inner arcs first, so an outer
     *    region always finds its inner arcs already scored.
     * 2. The top level is aligned inside an imagined arc spanning both
     *    sequences, i.e. structure boundaries (0, lenA+1) and (0, lenB+1).
     *
     * The optimal score ends up at M(lenA, lenB, lenA, lenB): both layers have
     * consumed both sequences completely.
     *
     * @return Optimal alignment score
     */
    score_t align();

    /**
     * @brief Reconstruct the alignments from the DP matrices
     *
     * Walks back from M(lenA, lenB, lenA, lenB), rebuilding the U and V layers,
     * the consensus structure, and the per-sequence dot-bracket strings. The
     * results are then available via get_alignment_U(), get_alignment_V(),
     * get_consensus_structure(), format_alignment(), and friends.
     *
     * @pre align() has been called.
     * @throws std::runtime_error if no predecessor reproduces the stored score,
     *         which would indicate an inconsistency between the forward pass and
     *         the traceback.
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

    std::string get_consensus_structure() const { return consensus_structure_; }
    std::string get_dot_bracket_A() const { return dot_bracket_A_; }
    std::string get_dot_bracket_B() const { return dot_bracket_B_; }

    /// Per-column shift annotation: '=' where U==V, 'S' where gap patterns differ.
    std::string get_shift_string() const {
        std::string result;
        size_t len = alignment_U_seqA_.size();
        result.reserve(len);
        for (size_t i = 0; i < len; ++i) {
            bool u_gap_a = (alignment_U_seqA_[i] == '-');
            bool u_gap_b = (alignment_U_seqB_[i] == '-');
            bool v_gap_a = (alignment_V_seqA_[i] == '-');
            bool v_gap_b = (alignment_V_seqB_[i] == '-');
            result += (u_gap_a == v_gap_a && u_gap_b == v_gap_b) ? '=' : 'S';
        }
        return result;
    }

    /// Consensus sequence: nucleotide where U_A and U_B agree, '.' otherwise.
    std::string get_consensus_sequence() const {
        std::string result;
        size_t len = alignment_U_seqA_.size();
        result.reserve(len);
        for (size_t i = 0; i < len; ++i) {
            char a = alignment_U_seqA_[i], b = alignment_U_seqB_[i];
            result += (a != '-' && b != '-' && a == b) ? a : '.';
        }
        return result;
    }

private:
    // Sequences
    std::unique_ptr<LocARNA::Sequence> seqA_;
    std::unique_ptr<LocARNA::Sequence> seqB_;

    // Preprocessing data
    std::unique_ptr<LocARNA::RnaData> rna_dataA_;
    std::unique_ptr<LocARNA::RnaData> rna_dataB_;
    std::unique_ptr<LocARNA::ArcMatches> arc_matches_;

    // Scoring
    // locarna_scoring_params_ MUST be declared before locarna_scoring_:
    // LocARNA::Scoring stores params by pointer, so the ScoringParams object
    // must outlive the Scoring object. Members are destroyed in reverse
    // declaration order, so declaring params first ensures it lives long enough. TODO: comment still needed?
    std::unique_ptr<LocARNA::ScoringParams> locarna_scoring_params_;
    std::unique_ptr<LocARNA::Scoring> locarna_scoring_;
    ShiftScoring shift_scoring_;

    // Heuristic parameters
    size_type max_shifts_;  ///< Maximum allowed number of shifts (δ_max)
    int verbose_;           ///< Verbosity level (0=silent, 1=basic, 2=detailed, 3=debug)

    // Dynamic programming matrices
    std::unique_ptr<ShiftMatrixM<score_t>> M_;  ///< Main DP matrix (4D)
    std::unique_ptr<ShiftMatrixD<score_t>> D_;  ///< Structure DP matrix (6D)

    // Results
    score_t alignment_score_;  ///< Computed optimal score

    // Alignment storage (populated by traceback)
    std::string alignment_U_seqA_;  ///< Sequence A in alignment U (with gaps)
    std::string alignment_U_seqB_;  ///< Sequence B in alignment U (with gaps)
    std::string alignment_V_seqA_;  ///< Sequence A in alignment V (with gaps)
    std::string alignment_V_seqB_;  ///< Sequence B in alignment V (with gaps)

    // Consensus structure annotation (populated by traceback)
    std::string consensus_structure_;  ///< Structure annotation: '(' and ')' for paired, '.' for unpaired
    std::vector<size_t> arc_stack_;    ///< Stack to track opening positions during traceback

    // Dot-bracket annotations aligned to structure layer (populated by traceback)
    std::string dot_bracket_A_;  ///< Dot-bracket centroid structure for seq A (prob>=0.5), aligned to V layer
    std::string dot_bracket_B_;  ///< Dot-bracket centroid structure for seq B (prob>=0.5), aligned to V layer

    // Per-position bracket character for each sequence (1-based, computed before traceback)
    std::vector<char> bracket_A_;  ///< bracket_A_[i] = '(' / ')' / '.' for position i in seqA
    std::vector<char> bracket_B_;  ///< bracket_B_[i] = '(' / ')' / '.' for position i in seqB

    // ========== Traceback helper functions (LocARNA-style recursive structure) ==========

    /**
     * @brief Core recursive traceback through M matrix
     *
     * Recursively traces back through the M matrix to reconstruct the optimal
     * bi-alignment. Tries all 15 unpaired column types (Case 1) and all arc
     * match cases (Case 2), finding which led to the current M value.
     *
     * @param al Left structure boundary for A
     * @param ar Right structure boundary for A
     * @param bl Left structure boundary for B
     * @param br Right structure boundary for B
     * @param x1 Sequence starting position for A
     * @param x2 Sequence starting position for B
     * @param y1 Current sequence position for A
     * @param y2 Current sequence position for B
     * @param y3 Current structure position for A
     * @param y4 Current structure position for B
     */
    void trace_M(size_type al, size_type ar,
                 size_type bl, size_type br,
                 size_type x1, size_type x2,
                 size_type y1, size_type y2,
                 size_type y3, size_type y4);

    /**
     * @brief Handle paired/arc match cases during traceback (Case 2)
     *
     * Called from trace_M when trying Case 2 (paired positions with arc match).
     * Iterates over all arc matches with right ends at (y3, y4) and tries to find
     * which one (along with which sequence offsets and gap patterns) led to current M value.
     *
     * @param al Left structure boundary for A
     * @param ar Right structure boundary for A
     * @param bl Left structure boundary for B
     * @param br Right structure boundary for B
     * @param x1 Sequence starting position for A
     * @param x2 Sequence starting position for B
     * @param y1 Current sequence position for A
     * @param y2 Current sequence position for B
     * @param y3 Current structure position for A
     * @param y4 Current structure position for B
     * @param current_score Score at M(y1, y2, y3, y4) for validation
     */
    void trace_M_paired(size_type al, size_type ar,
                       size_type bl, size_type br,
                       size_type x1, size_type x2,
                       size_type y1, size_type y2,
                       size_type y3, size_type y4,
                       score_t current_score);

    /**
     * @brief Traceback inside a matched arc pair
     *
     * Traces the region inside matched arcs. Since D = M(y1,y2,ar-1,br-1) + arcmatch,
     * this simply calls trace_M for the interior region.
     *
     * @param arcA Matched arc from sequence A
     * @param arcB Matched arc from sequence B
     * @param x1 Sequence position at left end of arcA
     * @param x2 Sequence position at left end of arcB
     * @param y1 Sequence position at right end of arcA
     * @param y2 Sequence position at right end of arcB
     */
    void trace_D(const LocARNA::BasePairs__Arc& arcA,
                 const LocARNA::BasePairs__Arc& arcB,
                 size_type x1, size_type x2,
                 size_type y1, size_type y2,
                 score_t arc_match_score);

    /**
     * @brief Helper: Append a column to the alignment
     *
     * Adds one column to the bi-alignment based on the column type.
     * Appends to alignment strings (forward construction).
     *
     * @param c1 1 if seqA advances in U layer, 0 for gap
     * @param c2 1 if seqB advances in U layer, 0 for gap
     * @param c3 1 if seqA advances in V layer, 0 for gap
     * @param c4 1 if seqB advances in V layer, 0 for gap
     * @param y1 Current sequence position in A
     * @param y2 Current sequence position in B
     * @param y3 Current structure position in A
     * @param y4 Current structure position in B
     * @param is_arc_left If true, this is the left position of an arc (opening bracket)
     * @param is_arc_right If true, this is the right position of an arc (closing bracket)
     */
    void append_column(int c1, int c2, int c3, int c4,
                      size_type y1, size_type y2,
                      size_type y3, size_type y4,
                      bool is_arc_left = false,
                      bool is_arc_right = false);

    // ========== Forward pass ==========

    /**
     * @brief Fill D for every arc match
     *
     * Iterates left arc endpoints (al, bl) in descending order, so inner arcs
     * are complete before any enclosing arc needs them. For each (al, bl) and
     * each valid sequence offset (x1, x2) within δ_max:
     * 1. align_in_arcmatch() fills a local M for the arc region,
     * 2. fill_D_entries() reads the D values out of it.
     *
     * @post D holds an entry for every arc match and endpoint-offset
     *       combination that is reachable within δ_max.
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
     * @param x1 Sequence position at left end in A (can shift from al+1)
     * @param x2 Sequence position at left end in B (can shift from bl+1)
     * @param max_y1 Maximum sequence position in A to fill M matrix up to
     * @param max_y2 Maximum sequence position in B to fill M matrix up to
     */
    void align_in_arcmatch(size_type al, size_type ar,
                           size_type bl, size_type br,
                           size_type x1, size_type x2,
                           size_type max_y1, size_type max_y2);

    /**
     * @brief Extract D matrix entries for arc matches with left ends (al, bl)
     *
     * For each arc match with left ends at (al, bl) and each valid
     * sequence offset combination:
     * D(am, z1, z2, y1, y2) = M(y1, y2, ar-1, br-1) + arcmatch_score(am) TODO
     *
     * @param al Left endpoint in sequence A (structure layer)
     * @param bl Left endpoint in sequence B (structure layer)
     * @param x1 Sequence position at left end in A
     * @param x2 Sequence position at left end in B
     */
    void fill_D_entries(size_type al, size_type bl,
                        size_type x1, size_type x2);



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
};

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFT_ALIGNER_HH
