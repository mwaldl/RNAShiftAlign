#include "shift_aligner.hh"

#include <sstream>
#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/pfold_params.hh>
#include <LocARNA/anchor_constraints.hh>
#include <LocARNA/trace_controller.hh>
#include <LocARNA/multiple_alignment.hh>

namespace RNAShiftAlign {

ShiftAligner::ShiftAligner(const std::string &seqA,
                           const std::string &seqB,
                           const PreprocessingParams &prep_params,
                           const ScoringParams &scoring_params)
    : shift_scoring_(score_t(scoring_params.delta)),
      max_shifts_(scoring_params.max_shifts),
      alignment_score_(score_t(0))
{
    // Step 1: Create sequences
    seqA_ = std::make_unique<LocARNA::Sequence>("seqA", seqA);
    seqB_ = std::make_unique<LocARNA::Sequence>("seqB", seqB);

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Step 2: Compute base pair probabilities via ViennaRNA partition function
    // Sequence → MultipleAlignment → RnaEnsemble → RnaData
    LocARNA::PFoldParams pfold_params;  // Use default ViennaRNA parameters

    // Create MultipleAlignment from single sequences
    auto ma_A = std::make_unique<LocARNA::MultipleAlignment>("seqA", seqA);
    auto ma_B = std::make_unique<LocARNA::MultipleAlignment>("seqB", seqB);

    // Create RnaEnsemble from MultipleAlignments
    auto rna_ensembleA = std::make_unique<LocARNA::RnaEnsemble>(
        *ma_A,
        pfold_params,
        false,  // no stacking
        false   // not in_loop
    );

    auto rna_ensembleB = std::make_unique<LocARNA::RnaEnsemble>(
        *ma_B,
        pfold_params,
        false,  // no stacking
        false   // not in_loop
    );

    // Create RnaData from ensembles
    rna_dataA_ = std::make_unique<LocARNA::RnaData>(
        *rna_ensembleA,
        prep_params.min_prob,
        prep_params.max_bps_length_ratio,
        pfold_params
    );

    rna_dataB_ = std::make_unique<LocARNA::RnaData>(
        *rna_ensembleB,
        prep_params.min_prob,
        prep_params.max_bps_length_ratio,
        pfold_params
    );

    // Step 3: Compute arc matches between the two sequences
    // Need to provide TraceController and AnchorConstraints (use defaults)

    // Default constraints: no anchors
    LocARNA::AnchorConstraints seq_constraints(lenA, "", lenB, "", true);

    // Default trace controller: no restrictions (use -1 for max_diff to disable, false for no relaxation)
    auto trace_controller = std::make_unique<LocARNA::TraceController>(
        *seqA_,
        *seqB_,
        nullptr,  // no reference alignment
        -1,       // max_diff = -1 (no restriction)
        false     // no relaxation
    );

    // Compute arc matches with proper parameters
    size_type max_diff_am_value = (prep_params.max_diff_am != -1)
        ? size_type(prep_params.max_diff_am)
        : std::max(lenA, lenB);
    size_type max_diff_at_am_value = (prep_params.max_diff_at_am != -1)
        ? size_type(prep_params.max_diff_at_am)
        : std::max(lenA, lenB);

    arc_matches_ = std::make_unique<LocARNA::ArcMatches>(
        *rna_dataA_,
        *rna_dataB_,
        prep_params.min_prob,
        max_diff_am_value,
        max_diff_at_am_value,
        *trace_controller,
        seq_constraints
    );

    // Step 4: Initialize scoring scheme using named argument pattern
    // Compute expected probabilities
    double exp_probA = rna_dataA_->arc_cutoff_prob();
    double exp_probB = rna_dataB_->arc_cutoff_prob();

    // Create scoring parameters using named argument pattern
    auto locarna_scoring_params = LocARNA::ScoringParams(
        LocARNA::ScoringParams::match(scoring_params.match),
        LocARNA::ScoringParams::mismatch(scoring_params.mismatch),
        LocARNA::ScoringParams::indel(scoring_params.indel),
        LocARNA::ScoringParams::indel_opening(scoring_params.indel_opening),
        LocARNA::ScoringParams::struct_weight(scoring_params.struct_weight),
        LocARNA::ScoringParams::tau_factor(scoring_params.tau_factor),
        LocARNA::ScoringParams::exp_probA(exp_probA),
        LocARNA::ScoringParams::exp_probB(exp_probB)
    );

    // Initialize LocARNA Scoring object (pass nullptr for match_probs - not needed for basic scoring)
    locarna_scoring_ = std::make_unique<LocARNA::Scoring>(
        *seqA_,
        *seqB_,
        *rna_dataA_,
        *rna_dataB_,
        *arc_matches_,
        nullptr,  // match_probs not needed for basic scoring
        locarna_scoring_params
    );

    // Step 5: Allocate shift-aware DP matrices
    // M matrix is allocated fresh in each align_in_arcmatch() call
    // since entries are independent between different starting columns

    // D matrix: 6D matrix for structure alignment D(arc_a, arc_b, x1, x2, y1, y2)
    // Dimensions: number of base pairs in each RNA
    size_type num_bps_A = arc_matches_->get_base_pairsA().num_bps();
    size_type num_bps_B = arc_matches_->get_base_pairsB().num_bps();
    D_ = std::make_unique<LocARNA::ShiftMatrixD<score_t>>(
        num_bps_A,
        num_bps_B,
        max_shifts_
    );
    // Note: Individual 4D offset matrices are created on-demand in fill_D()
}

ShiftAligner::score_t
ShiftAligner::align() {
    // Two-phase Sankoff algorithm following LocARNA pattern:
    // 1. Fill D-matrix (processes arcs in descending order of left endpoints)
    // 2. Align top level with imagined surrounding arc (0, lenA+1, 0, lenB+1)

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Phase 1: Fill D matrix for all arc matches
    align_D();

    // Phase 2: Top-level alignment with imagined arc around entire sequences
    // Structure boundaries: (0, lenA+1, 0, lenB+1) - virtual arc endpoints
    // Sequence boundaries: (0, 0) to (lenA, lenB) - all positions used
    align_in_arcmatch(0, lenA + 1, 0, lenB + 1,
                      0, 0,           // x1, x2: sequence start
                      lenA, lenB);    // y1, y2: sequence end

    // Optimal score is at M(lenA, lenB, lenA, lenB)
    // Both sequence and structure layers must consume all positions
    alignment_score_ = M_->get(lenA, lenB, lenA, lenB);
    return alignment_score_;
}

void
ShiftAligner::traceback() {
    // Reconstruct optimal bi-alignment from M matrix
    // Start from endpoint and work backwards to (0,0,0,0)

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Clear any previous alignments
    alignment_U_seqA_.clear();
    alignment_U_seqB_.clear();
    alignment_V_seqA_.clear();
    alignment_V_seqB_.clear();

    // Define column types (same as in fill_M_unpaired)
    struct Column {
        int c1, c2, c3, c4;
    };

    // Use same column ordering as forward pass to ensure consistency
    // No-shift columns first, then shift columns
    const std::vector<Column> valid_columns = {
        // No-shift columns (3 total): c_U == c_V
        {1,1,1,1},  // match-match
        {1,0,1,0},  // del_A-del_A
        {0,1,0,1},  // ins_A-ins_A

        // Shift columns (12 total): c_U != c_V
        {1,1,1,0}, {1,1,0,1}, {1,1,0,0},
        {1,0,1,1}, {1,0,0,1}, {1,0,0,0},
        {0,1,1,1}, {0,1,1,0}, {0,1,0,0},
        {0,0,1,1}, {0,0,1,0}, {0,0,0,1}
    };

    // Start from endpoint
    size_type y1 = lenA;
    size_type y2 = lenB;
    size_type y3 = lenA;
    size_type y4 = lenB;

    // Traceback until we reach (0,0,0,0)
    while (y1 > 0 || y2 > 0 || y3 > 0 || y4 > 0) {
        score_t current_score = M_->get(y1, y2, y3, y4);

        // Find which column type led to this cell
        bool found = false;
        for (const auto& col : valid_columns) {
            // Compute predecessor
            int x1 = y1 - col.c1;
            int x2 = y2 - col.c2;
            int x3 = y3 - col.c3;
            int x4 = y4 - col.c4;

            // Check validity
            if (x1 < 0 || x2 < 0 || x3 < 0 || x4 < 0)
                continue;

            // Check δ_max constraint
            if (std::abs(x1 - x3) > static_cast<int>(max_shifts_))
                continue;
            if (std::abs(x2 - x4) > static_cast<int>(max_shifts_))
                continue;

            // Get predecessor score
            score_t pred_score = M_->get(x1, x2, x3, x4);

            // Compute column score (same logic as forward)
            score_t u_s(0);
            ColumnType col_U;
            if (col.c1 == 1 && col.c2 == 1) {
                u_s = score_t(locarna_scoring_->basematch(y1, y2));
                col_U = ColumnType::MATCH;
            } else if (col.c1 == 1 && col.c2 == 0) {
                u_s = score_t(locarna_scoring_->gapA(y1));
                col_U = ColumnType::DEL_A;
            } else if (col.c1 == 0 && col.c2 == 1) {
                u_s = score_t(locarna_scoring_->gapB(y2));
                col_U = ColumnType::INS_A;
            } else {
                col_U = ColumnType::MATCH;
            }

            score_t v_s(0);
            ColumnType col_V;
            if (col.c3 == 1 && col.c4 == 1) {
                v_s = score_t(locarna_scoring_->basematch(y3, y4));
                col_V = ColumnType::MATCH;
            } else if (col.c3 == 1 && col.c4 == 0) {
                v_s = score_t(locarna_scoring_->gapA(y3));
                col_V = ColumnType::DEL_A;
            } else if (col.c3 == 0 && col.c4 == 1) {
                v_s = score_t(locarna_scoring_->gapB(y4));
                col_V = ColumnType::INS_A;
            } else {
                col_V = ColumnType::MATCH;
            }

            score_t w_s = shift_scoring_.shift_penalty(col_U, col_V);
            score_t column_score = u_s + v_s + w_s;

            // Check if this is the path we took
            if (pred_score + column_score == current_score) {
                // Found the traceback step!
                // Append to alignment strings (we're building backwards, so prepend)

                // For U layer (sequence alignment)
                // Note: LocARNA sequences are 1-indexed (position 0 contains a space)
                if (col.c1 == 1) {
                    alignment_U_seqA_ = seqA_->seqentry(0).seq()[y1] + alignment_U_seqA_;
                } else {
                    alignment_U_seqA_ = "-" + alignment_U_seqA_;
                }

                if (col.c2 == 1) {
                    alignment_U_seqB_ = seqB_->seqentry(0).seq()[y2] + alignment_U_seqB_;
                } else {
                    alignment_U_seqB_ = "-" + alignment_U_seqB_;
                }

                // For V layer (structure alignment)
                if (col.c3 == 1) {
                    alignment_V_seqA_ = seqA_->seqentry(0).seq()[y3] + alignment_V_seqA_;
                } else {
                    alignment_V_seqA_ = "-" + alignment_V_seqA_;
                }

                if (col.c4 == 1) {
                    alignment_V_seqB_ = seqB_->seqentry(0).seq()[y4] + alignment_V_seqB_;
                } else {
                    alignment_V_seqB_ = "-" + alignment_V_seqB_;
                }

                // Move to predecessor
                y1 = x1;
                y2 = x2;
                y3 = x3;
                y4 = x4;

                found = true;
                break;
            }
        }

        if (!found) {
            throw std::runtime_error("Traceback failed: no valid predecessor found");
        }
    }
}

std::string
ShiftAligner::format_alignment() const {
    std::ostringstream oss;

    // Helper to create structure annotation (currently all dots for unpaired)
    auto make_structure = [](const std::string& seq) {
        std::string structure;
        for (char c : seq) {
            structure += (c == '-') ? '-' : '.';
        }
        return structure;
    };

    // Helper to identify shift positions
    auto find_shifts = [](const std::string& u_A, const std::string& u_B,
                          const std::string& v_A, const std::string& v_B) {
        std::string shifts;
        size_t len = u_A.length();  // All should be same length
        for (size_t i = 0; i < len; ++i) {
            // Compare column patterns: (u_A[i], u_B[i]) vs (v_A[i], v_B[i])
            bool u_both = (u_A[i] != '-' && u_B[i] != '-');
            bool u_gap_A = (u_A[i] != '-' && u_B[i] == '-');
            bool u_gap_B = (u_A[i] == '-' && u_B[i] != '-');
            bool u_both_gap = (u_A[i] == '-' && u_B[i] == '-');

            bool v_both = (v_A[i] != '-' && v_B[i] != '-');
            bool v_gap_A = (v_A[i] != '-' && v_B[i] == '-');
            bool v_gap_B = (v_A[i] == '-' && v_B[i] != '-');
            bool v_both_gap = (v_A[i] == '-' && v_B[i] == '-');

            // Shift occurs when gap patterns differ
            bool is_shift = !(u_both == v_both && u_gap_A == v_gap_A &&
                             u_gap_B == v_gap_B && u_both_gap == v_both_gap);

            shifts += is_shift ? 'S' : '=';
        }
        return shifts;
    };

    // Create structure annotations
    std::string strA_U = make_structure(alignment_U_seqA_);
    std::string strB_U = make_structure(alignment_U_seqB_);
    std::string strA_V = make_structure(alignment_V_seqA_);
    std::string strB_V = make_structure(alignment_V_seqB_);

    // Find shift positions
    std::string shifts = find_shifts(alignment_U_seqA_, alignment_U_seqB_,
                                     alignment_V_seqA_, alignment_V_seqB_);

    // Format output
    oss << "Sequence Alignment (U):\n";
    oss << "strA: " << strA_U << "\n";
    oss << "seqA: " << alignment_U_seqA_ << "\n";
    oss << "seqB: " << alignment_U_seqB_ << "\n";
    oss << "strB: " << strB_U << "\n";
    oss << "\n";

    oss << "Structure Alignment (V):\n";
    oss << "seqA: " << alignment_V_seqA_ << "\n";
    oss << "strA: " << strA_V << "\n";
    oss << "strB: " << strB_V << "\n";
    oss << "seqB: " << alignment_V_seqB_ << "\n";
    oss << "\n";

    oss << "Shift Annotation:\n";
    oss << "      " << shifts << "\n";
    oss << "      (= indicates U==V, S indicates shift)\n";

    return oss.str();
}





// ============================================================================
// NEW IMPLEMENTATION: LocARNA-style structure following program_flow_analysis.md
// ============================================================================
//
// These functions implement the restructured algorithm that properly handles
// shifts between sequence and structure layers. The key insight is that
// arc endpoints have fixed structure positions (al, bl, ar, br) but the
// corresponding sequence positions (x1, x2, y1, y2) can shift within δ_max.
//
// Program flow:
// 1. init_D() - Create D matrix storage
// 2. align_D() - Fill D matrix for all arc matches
//    - For each (al, bl) in descending order:
//      - For each valid (x1, x2) offset:
//        - align_in_arcmatch() to fill local M
//        - fill_D_entries() to extract D values
// 3. align() calls align_in_arcmatch for top level (imagined surrounding arc)
//
// See planning_docs/program_flow_analysis.md for full details.
// ============================================================================

void
ShiftAligner::init_D() {
    // ------------------------------------------------------------
    // Initialize D matrix storage
    //
    // D matrix allocation is done in the constructor:
    // - D_ is created with dimensions (num_bps_A, num_bps_B, max_shifts_)
    //
    // This function can be used to initialize D matrix values if needed.
    // Currently, D values are filled on-demand by fill_D_entries().
    //
    // The offset matrices for each arc pair are created when needed
    // via D_->create_offsetmatrix().
    // ------------------------------------------------------------

    // D matrix is already allocated in constructor
    // Individual offset matrices are created on-demand in fill_D_entries()
}

void
ShiftAligner::align_D() {
    // ------------------------------------------------------------
    // Fill D matrix for all arc matches
    //
    // Following LocARNA's align_D() pattern but with sequence offsets.
    //
    // Outer loop: iterate over structure positions (al, bl) in descending order
    //   This ensures inner arcs are processed before outer arcs.
    //
    // Middle loop: iterate over valid sequence offsets (x1, x2)
    //   x1 can range from al - max_shifts to al + max_shifts
    //   x2 can range from bl - max_shifts to bl + max_shifts
    //
    // For each combination:
    //   1. Find max right ends for arcs starting at (al, bl)
    //   2. Compute corresponding y1, y2 sequence positions
    //   3. Call align_in_arcmatch(al, ar, bl, br, x1, x2, y1, y2)
    //   4. Call fill_D_entries(al, bl, x1, x2)
    //
    // Key difference from LocARNA:
    // - We must iterate over x1, x2 offsets at each (al, bl)
    // - The y1, y2 offsets are computed based on the arc lengths
    //   and the starting offsets
    //
    // Implementation steps:
    // 1. for al = lenA ... 1 (descending)
    // 2.   for bl = lenB ... 1 (descending)
    // 3.     Get arcmatches at (al, bl) using common_left_end_list
    // 4.     If empty, continue
    // 5.     Find max_ar, max_br for arcs starting here
    // 6.     for shift_x1 = -max_shifts ... +max_shifts
    // 7.       for shift_x2 = -max_shifts ... +max_shifts
    // 8.         x1 = al + shift_x1, x2 = bl + shift_x2
    // 9.         Check bounds
    // 10.        Compute y1, y2 from ar, br and same relative shifts
    // 11.        align_in_arcmatch(al, max_ar, bl, max_br, x1, x2, y1, y2)
    // 12.        fill_D_entries(al, bl, x1, x2)
    // ------------------------------------------------------------

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Iterate over left endpoints in descending order
    for (size_type al = lenA; al > 0; --al) {
        for (size_type bl = lenB; bl > 0; --bl) {

            // Get arc matches with left ends at (al, bl)
            const auto& arcmatches_at_left = arc_matches_->common_left_end_list(al, bl);
            if (arcmatches_at_left.empty())
                continue;

            // Find maximal right ends
            size_type max_ar = al;
            size_type max_br = bl;
            for (const auto& am_idx : arcmatches_at_left) {
                const auto& am = arc_matches_->arcmatch(am_idx);
                max_ar = std::max(max_ar, am.arcA().right());
                max_br = std::max(max_br, am.arcB().right());
            }

            // Skip if no valid arcs
            if (max_ar == al || max_br == bl)
                continue;

            // Iterate over sequence offsets at left end (x1, x2)
            for (int shift_x1 = -static_cast<int>(max_shifts_);
                 shift_x1 <= static_cast<int>(max_shifts_); ++shift_x1) {
                for (int shift_x2 = -static_cast<int>(max_shifts_);
                     shift_x2 <= static_cast<int>(max_shifts_); ++shift_x2) {

                    int x1 = static_cast<int>(al) + shift_x1;
                    int x2 = static_cast<int>(bl) + shift_x2;

                    // Bounds check for x1, x2
                    if (x1 < 1 || x2 < 1)
                        continue;
                    if (x1 > static_cast<int>(lenA) || x2 > static_cast<int>(lenB))
                        continue;

                    // Compute maximum y1, y2 (right end sequence positions)
                    // Fill M up to max possible shift from structure right ends
                    size_type max_y1 = std::min(max_ar + max_shifts_, lenA);
                    size_type max_y2 = std::min(max_br + max_shifts_, lenB);

                    // Fill M matrix for this arc region with all possible right ends
                    align_in_arcmatch(al, max_ar, bl, max_br,
                                      static_cast<size_type>(x1),
                                      static_cast<size_type>(x2),
                                      max_y1, max_y2);

                    // Extract D entries for all valid (y1, y2) combinations
                    fill_D_entries(al, bl,
                                   static_cast<size_type>(x1),
                                   static_cast<size_type>(x2));
                }
            }
        }
    }
}

void
ShiftAligner::align_in_arcmatch(size_type al, size_type ar,
                                 size_type bl, size_type br,
                                 size_type x1, size_type x2,
                                 size_type y1, size_type y2) {
    // ------------------------------------------------------------
    // Fill M matrix for region inside an arc match
    //
    // This corresponds to LocARNA's align_in_arcmatch but handles
    // both sequence and structure positions.
    //
    // Parameters:
    // - (al, ar, bl, br): Structure layer boundaries (fixed by arc)
    // - (x1, x2): Sequence positions at left end (can shift from al, bl)
    // - (y1, y2): Sequence positions at right end (can shift from ar, br)
    //
    // IMPORTANT: Each call creates a fresh M matrix for this specific
    // (x1, x2, al, bl) starting column. The M matrix entries are
    // independent between different starting columns.
    //
    // Steps:
    // 1. Create fresh M matrix and fill with -infinity
    // 2. Initialize M matrix boundaries: init_M(...)
    // 3. For each structure position (y3, y4) in [al+1, ar-1] x [bl+1, br-1]:
    //    For each valid sequence position (seq_y1, seq_y2):
    //      M(seq_y1, seq_y2, y3, y4) = align_noex(...)
    //
    // Future optimization: Only allocate the slice between
    // (x1, x2, al, bl) and (y1, y2, ar, br) instead of the full matrix.
    // ------------------------------------------------------------

    // Create fresh M matrix for this arc region
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();
    M_ = std::make_unique<LocARNA::ShiftMatrixM<score_t>>(
        lenA + 1,
        lenB + 1,
        max_shifts_
    );
    M_->fill(LocARNA::infty_score_t::neg_infty);
    

    // Initialize boundaries
    init_M(al, ar, bl, br, x1, x2, y1, y2);

    // Fill M matrix for interior positions
    for (size_type y3 = al + 1; y3 < ar; ++y3) {
        for (size_type y4 = bl + 1; y4 < br; ++y4) {

            // For each valid sequence position around (y3, y4)
            for (int shift1 = -static_cast<int>(max_shifts_);
                 shift1 <= static_cast<int>(max_shifts_); ++shift1) {
                for (int shift2 = -static_cast<int>(max_shifts_);
                     shift2 <= static_cast<int>(max_shifts_); ++shift2) {

                    int seq_y1 = static_cast<int>(y3) + shift1;
                    int seq_y2 = static_cast<int>(y4) + shift2;

                    // Bounds check
                    if (seq_y1 < static_cast<int>(x1) || seq_y1 > static_cast<int>(y1))
                        continue;
                    if (seq_y2 < static_cast<int>(x2) || seq_y2 > static_cast<int>(y2))
                        continue;

                    // Compute M value using core recursion
                    score_t score = align_noex(al, bl, x1, x2,
                                               static_cast<size_type>(seq_y1),
                                               static_cast<size_type>(seq_y2),
                                               y3, y4);

                    M_->set(static_cast<size_type>(seq_y1),
                            static_cast<size_type>(seq_y2),
                            y3, y4, score);
                }
            }
        }
    }
}

ShiftAligner::score_t
ShiftAligner::align_noex(size_type al, size_type bl,
                          size_type x1, size_type x2,
                          size_type y1, size_type y2,
                          size_type y3, size_type y4) {
    // ------------------------------------------------------------
    // Core recursion: compute optimal M value at position
    //
    // This is the heart of the algorithm. It optimizes over:
    // - Case 1: All 15 unpaired column types
    // - Case 2: All arc matches with right ends at (y3, y4)
    //
    // Parameters:
    // - (al, bl): Left boundaries of current arc region (structure)
    // - (x1, x2): Starting sequence positions
    // - (y1, y2): Current sequence position
    // - (y3, y4): Current structure position
    //
    // Returns: Optimal score for M(y1, y2, y3, y4)
    //
    // Implementation:
    // 1. Initialize max_score = -infinity
    // 2. Case 1: max_score = max(max_score, compute_unpaired_score(...))
    // 3. Case 2: max_score = max(max_score, compute_paired_score(...))
    // 4. Return max_score
    // ------------------------------------------------------------

    score_t max_score = score_t::neg_infty;

    // Case 1: Unpaired positions
    max_score = std::max(max_score, compute_unpaired_score(al, bl, x1, x2, y1, y2, y3, y4));

    // Case 2: Paired positions (arc matches)
    max_score = std::max(max_score, compute_paired_score(al, bl, x1, x2, y1, y2, y3, y4));

    return max_score;
}

void
ShiftAligner::fill_D_entries(size_type al, size_type bl,
                              size_type x1, size_type x2) {
    // ------------------------------------------------------------
    // Extract D matrix entries for arc matches with left ends (al, bl)
    //
    // For each arc match with left ends at (al, bl):
    //   ar = arcA.right(), br = arcB.right()
    //
    //   For each valid (y1, y2) sequence position at right end:
    //     D(am, x1, x2, y1, y2) = M(y1, y2, ar-1, br-1) + arcmatch_score(am)
    //
    // Note: x1, x2 are the sequence positions at the left end (passed in)
    //       y1, y2 are iterated over based on ar, br and max_shifts
    //
    // The M matrix has been filled up to the maximum possible y1, y2
    // so all valid entries are available for extraction.
    // ------------------------------------------------------------

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    const auto& arcmatches_at_left = arc_matches_->common_left_end_list(al, bl);

    for (const auto& am_idx : arcmatches_at_left) {
        const auto& am = arc_matches_->arcmatch(am_idx);
        const auto& arcA = am.arcA();
        const auto& arcB = am.arcB();

        size_type ar = arcA.right();
        size_type br = arcB.right();

        // Create offset matrix for this arc pair if needed
        D_->create_offsetmatrix(arcA.idx(), arcB.idx());

        // Structure positions just inside arc right ends
        size_type y3_fixed = ar - 1;
        size_type y4_fixed = br - 1;

        // Iterate over all valid sequence positions at right end
        // y1, y2 can shift within δ_max of the structure positions
        for (int shift_y1 = -static_cast<int>(max_shifts_);
             shift_y1 <= static_cast<int>(max_shifts_); ++shift_y1) {
            for (int shift_y2 = -static_cast<int>(max_shifts_);
                 shift_y2 <= static_cast<int>(max_shifts_); ++shift_y2) {

                int y1 = static_cast<int>(y3_fixed) + shift_y1;
                int y2 = static_cast<int>(y4_fixed) + shift_y2;

                // Bounds check - y1, y2 must be inside the arc region
                // and within sequence bounds
                if (y1 <= static_cast<int>(al) || y2 <= static_cast<int>(bl))
                    continue;
                if (y1 > static_cast<int>(lenA) || y2 > static_cast<int>(lenB))
                    continue;
                // Also check they don't exceed the arc right ends
                // (sequence can extend past structure but not needed here)
                if (y1 >= static_cast<int>(ar) + static_cast<int>(max_shifts_) ||
                    y2 >= static_cast<int>(br) + static_cast<int>(max_shifts_))
                    continue;

                // Get M value - this is the score from (x1, x2, al, bl) to (y1, y2, y3_fixed, y4_fixed)
                score_t m_score = M_->get(static_cast<size_type>(y1),
                                          static_cast<size_type>(y2),
                                          y3_fixed, y4_fixed);

                // Skip if M value was never computed (still -infinity)
                if (m_score == score_t::neg_infty)
                    continue;

                // Compute D value = M + arcmatch_score
                score_t arc_score = score_t(locarna_scoring_->arcmatch(am));
                score_t d_score = m_score + arc_score;

                // Store in D matrix indexed by (arcA, arcB, x1, x2, y1, y2)
                D_->set(arcA, arcB, x1, x2,
                        static_cast<size_type>(y1),
                        static_cast<size_type>(y2),
                        d_score);
            }
        }
    }
}

void
ShiftAligner::init_M(size_type al, size_type ar,
                      size_type bl, size_type br,
                      size_type x1, size_type x2,
                      size_type y1, size_type y2) {
    // ------------------------------------------------------------
    // Initialize M matrix boundaries for arc region
    //
    // Sets up M matrix entries at the boundaries of the region.
    // This corresponds to LocARNA's init_state function.
    //
    // Key insight for 4D matrix initialization:
    // - Base case: M(x1, x2, al, bl) = 0
    // - Boundaries involve gaps in one or both layers
    //
    // We need to initialize all positions where at least one index
    // is at the boundary (al or bl for structure, x1 or x2 for sequence).
    // ------------------------------------------------------------

    // Base case: starting position
    M_->set(x1, x2, al, bl, score_t(0));

    // Initialize boundaries where structure B is at left edge (y4 = bl)
    // This means we're building up gaps in sequence B
    for (size_type y3 = al; y3 < ar; ++y3) {
        // For each valid sequence position y1 around y3
        for (int shift1 = -static_cast<int>(max_shifts_);
             shift1 <= static_cast<int>(max_shifts_); ++shift1) {

            int seq_y1 = static_cast<int>(y3) + shift1;

            // Bounds check
            if (seq_y1 < static_cast<int>(x1) || seq_y1 > static_cast<int>(y1))
                continue;

            // At boundary y4=bl, y2=x2: only gaps in B accumulated
            // Score = gap costs for A from x1 to seq_y1 + shift penalties
            //
            // The gap pattern is: sequence A advances, structure A advances,
            // sequence B stays at x2, structure B stays at bl

            if (seq_y1 == static_cast<int>(x1) && y3 == al) {
                // Base case already set
                continue;
            }

            // Compute accumulated gap cost for first row
            // M(seq_y1, x2, y3, bl) = M(seq_y1-1, x2, y3-1, bl) + gap_A + gap_A + shift
            // or handle the 15 column types that reach this position

            score_t row_score = score_t::neg_infty;

            // Try column types that keep y2=x2 and y4=bl
            // Valid types: (1,0,1,0) - gap in B for both layers
            //              (1,0,0,0) - gap in B for seq, gap in both for struct
            //              (0,0,1,0) - gap in both for seq, gap in B for struct

            // Type (1,0,1,0): advance sequence A and structure A
            if (seq_y1 > static_cast<int>(x1) && y3 > al) {
                int prev_y1 = seq_y1 - 1;
                size_type prev_y3 = y3 - 1;

                if (valid_shift(prev_y1, prev_y3)) {
                    score_t pred = M_->get(prev_y1, x2, prev_y3, bl);
                    score_t u_s = score_t(locarna_scoring_->gapA(seq_y1));
                    score_t v_s = score_t(locarna_scoring_->gapA(y3));
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::DEL_A, ColumnType::DEL_A);
                    row_score = std::max(row_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            // Type (1,0,0,0): advance sequence A only
            if (seq_y1 > static_cast<int>(x1)) {
                int prev_y1 = seq_y1 - 1;

                if (valid_shift(prev_y1, y3)) {
                    score_t pred = M_->get(prev_y1, x2, y3, bl);
                    score_t u_s = score_t(locarna_scoring_->gapA(seq_y1));
                    score_t v_s = score_t(0);  // structure gap in both
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::DEL_A, ColumnType::MATCH);
                    row_score = std::max(row_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            // Type (0,0,1,0): advance structure A only
            if (y3 > al) {
                size_type prev_y3 = y3 - 1;

                if (valid_shift(seq_y1, prev_y3)) {
                    score_t pred = M_->get(seq_y1, x2, prev_y3, bl);
                    score_t u_s = score_t(0);  // sequence gap in both
                    score_t v_s = score_t(locarna_scoring_->gapA(y3));
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A);
                    row_score = std::max(row_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            if (row_score > score_t::neg_infty) {
                M_->set(seq_y1, x2, y3, bl, row_score);
            }
        }
    }

    // Initialize boundaries where structure A is at left edge (y3 = al)
    // This means we're building up gaps in sequence A
    for (size_type y4 = bl; y4 < br; ++y4) {
        // For each valid sequence position y2 around y4
        for (int shift2 = -static_cast<int>(max_shifts_);
             shift2 <= static_cast<int>(max_shifts_); ++shift2) {

            int seq_y2 = static_cast<int>(y4) + shift2;

            // Bounds check
            if (seq_y2 < static_cast<int>(x2) || seq_y2 > static_cast<int>(y2))
                continue;

            if (seq_y2 == static_cast<int>(x2) && y4 == bl) {
                // Base case already set
                continue;
            }

            score_t col_score = score_t::neg_infty;

            // Try column types that keep y1=x1 and y3=al
            // Valid types: (0,1,0,1) - gap in A for both layers
            //              (0,1,0,0) - gap in A for seq, gap in both for struct
            //              (0,0,0,1) - gap in both for seq, gap in A for struct

            // Type (0,1,0,1): advance sequence B and structure B
            if (seq_y2 > static_cast<int>(x2) && y4 > bl) {
                int prev_y2 = seq_y2 - 1;
                size_type prev_y4 = y4 - 1;

                if (valid_shift(prev_y2, prev_y4)) {
                    score_t pred = M_->get(x1, prev_y2, al, prev_y4);
                    score_t u_s = score_t(locarna_scoring_->gapB(seq_y2));
                    score_t v_s = score_t(locarna_scoring_->gapB(y4));
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::INS_A, ColumnType::INS_A);
                    col_score = std::max(col_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            // Type (0,1,0,0): advance sequence B only
            if (seq_y2 > static_cast<int>(x2)) {
                int prev_y2 = seq_y2 - 1;

                if (valid_shift(prev_y2, y4)) {
                    score_t pred = M_->get(x1, prev_y2, al, y4);
                    score_t u_s = score_t(locarna_scoring_->gapB(seq_y2));
                    score_t v_s = score_t(0);  // structure gap in both
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::INS_A, ColumnType::MATCH);
                    col_score = std::max(col_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            // Type (0,0,0,1): advance structure B only
            if (y4 > bl) {
                size_type prev_y4 = y4 - 1;

                if (valid_shift(seq_y2, prev_y4)) {
                    score_t pred = M_->get(x1, seq_y2, al, prev_y4);
                    score_t u_s = score_t(0);  // sequence gap in both
                    score_t v_s = score_t(locarna_scoring_->gapB(y4));
                    score_t w_s = shift_scoring_.shift_penalty(ColumnType::MATCH, ColumnType::INS_A);
                    col_score = std::max(col_score, score_t(pred + u_s + v_s + w_s));
                }
            }

            if (col_score > score_t::neg_infty) {
                M_->set(x1, seq_y2, al, y4, col_score);
            }
        }
    }
}

ShiftAligner::score_t
ShiftAligner::compute_unpaired_score(size_type al, size_type bl,
                                      size_type x1, size_type x2,
                                      size_type y1, size_type y2,
                                      size_type y3, size_type y4) {
    // ------------------------------------------------------------
    // Compute best score for unpaired columns (Case 1)
    //
    // Iterates over all 15 valid column types and returns the best score.
    // Extracted from fill_M_local Case 1 logic.
    // ------------------------------------------------------------

    // All 15 valid column types (excluding (0,0,0,0))
    // Ordered to prefer no-shift columns first
    struct Column {
        int c1, c2, c3, c4;
    };

    const std::vector<Column> valid_columns = {
        // No-shift columns (3 total): c_U == c_V
        {1,1,1,1},  // match-match
        {1,0,1,0},  // del_A-del_A
        {0,1,0,1},  // ins_A-ins_A

        // Shift columns (12 total): c_U != c_V
        {1,1,1,0}, {1,1,0,1}, {1,1,0,0},
        {1,0,1,1}, {1,0,0,1}, {1,0,0,0},
        {0,1,1,1}, {0,1,1,0}, {0,1,0,0},
        {0,0,1,1}, {0,0,1,0}, {0,0,0,1}
    };

    score_t best_score = score_t::neg_infty;

    for (const auto& col : valid_columns) {
        // Compute predecessor position x = y - c
        int px1 = static_cast<int>(y1) - col.c1;
        int px2 = static_cast<int>(y2) - col.c2;
        int px3 = static_cast<int>(y3) - col.c3;
        int px4 = static_cast<int>(y4) - col.c4;

        // Check if predecessor is within region
        if (px1 < static_cast<int>(x1) || px2 < static_cast<int>(x2) ||
            px3 < static_cast<int>(al) || px4 < static_cast<int>(bl))
            continue;

        // Check δ_max constraint at predecessor
        if (std::abs(px1 - px3) > static_cast<int>(max_shifts_))
            continue;
        if (std::abs(px2 - px4) > static_cast<int>(max_shifts_))
            continue;

        // Get score from predecessor
        score_t pred_score = M_->get(px1, px2, px3, px4);

        // Compute column score s(y, c) = u_s + v_s + w_s

        // u_s: sequence alignment score for U layer (positions y1, y2)
        score_t u_s(0);
        ColumnType col_U;
        if (col.c1 == 1 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->basematch(y1, y2));
            col_U = ColumnType::MATCH;
        } else if (col.c1 == 1 && col.c2 == 0) {
            u_s = score_t(locarna_scoring_->gapA(y1));
            col_U = ColumnType::DEL_A;
        } else if (col.c1 == 0 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->gapB(y2));
            col_U = ColumnType::INS_A;
        } else {
            col_U = ColumnType::MATCH;  // (0,0) case
        }

        // v_s: structure score for V layer (positions y3, y4)
        score_t v_s(0);
        ColumnType col_V;
        if (col.c3 == 1 && col.c4 == 1) {
            v_s = score_t(locarna_scoring_->basematch(y3, y4));
            col_V = ColumnType::MATCH;
        } else if (col.c3 == 1 && col.c4 == 0) {
            v_s = score_t(locarna_scoring_->gapA(y3));
            col_V = ColumnType::DEL_A;
        } else if (col.c3 == 0 && col.c4 == 1) {
            v_s = score_t(locarna_scoring_->gapB(y4));
            col_V = ColumnType::INS_A;
        } else {
            col_V = ColumnType::MATCH;  // (0,0) case
        }

        // w_s: shift penalty
        score_t w_s = shift_scoring_.shift_penalty(col_U, col_V);

        score_t total_score = pred_score + u_s + v_s + w_s;

        if (total_score > best_score) {
            best_score = total_score;
        }
    }

    return best_score;
}

ShiftAligner::score_t
ShiftAligner::compute_paired_score(size_type al, size_type bl,
                                    size_type x1, size_type x2,
                                    size_type y1, size_type y2,
                                    size_type y3, size_type y4) {
    // ------------------------------------------------------------
    // Compute best score for paired columns (Case 2)
    //
    // Iterates over all arc matches with right ends at (y3, y4).
    // Extracted from fill_M_local Case 2 logic.
    // ------------------------------------------------------------

    score_t best_score = score_t::neg_infty;

    // Get all arc matches with right ends at current structure positions (y3, y4)
    const auto& arcmatches_at_right = arc_matches_->common_right_end_list(y3, y4);

    for (const auto& am_idx : arcmatches_at_right) {
        const auto& arc_match = arc_matches_->arcmatch(am_idx);
        const auto& arcA = arc_match.arcA();
        const auto& arcB = arc_match.arcB();

        size_type z3 = arcA.left();   // Left arc endpoint in structure A
        size_type z4 = arcB.left();   // Left arc endpoint in structure B

        // Arc must be inside current region
        if (z3 <= al || z4 <= bl)
            continue;

        // Iterate over all possible sequence positions z1, z2 at arc left ends
        for (int shift_z1 = -static_cast<int>(max_shifts_);
             shift_z1 <= static_cast<int>(max_shifts_); ++shift_z1) {
            for (int shift_z2 = -static_cast<int>(max_shifts_);
                 shift_z2 <= static_cast<int>(max_shifts_); ++shift_z2) {

                int z1 = static_cast<int>(z3) + shift_z1;
                int z2 = static_cast<int>(z4) + shift_z2;

                // Check bounds for z positions
                if (z1 < static_cast<int>(x1) || z2 < static_cast<int>(x2))
                    continue;

                // Iterate over valid gap patterns for z-column (arc left ends)
                for (int c1_z = 0; c1_z <= 1; ++c1_z) {
                    for (int c2_z = 0; c2_z <= 1; ++c2_z) {

                        // Iterate over valid gap patterns for y-column (arc right ends)
                        for (int c1_y = 0; c1_y <= 1; ++c1_y) {
                            for (int c2_y = 0; c2_y <= 1; ++c2_y) {

                                // M score before arc left endpoints
                                int px1 = z1 - c1_z;
                                int px2 = z2 - c2_z;
                                int px3 = static_cast<int>(z3) - 1;
                                int px4 = static_cast<int>(z4) - 1;

                                // Check bounds
                                if (px1 < static_cast<int>(x1) || px2 < static_cast<int>(x2) ||
                                    px3 < static_cast<int>(al) || px4 < static_cast<int>(bl))
                                    continue;

                                // Check δ_max constraint at predecessor
                                if (std::abs(px1 - px3) > static_cast<int>(max_shifts_))
                                    continue;
                                if (std::abs(px2 - px4) > static_cast<int>(max_shifts_))
                                    continue;

                                // Get M score before the arc
                                score_t m_before = M_->get(px1, px2, px3, px4);

                                // D score: optimal alignment inside matched arcs
                                score_t d_score = D_->get(arcA, arcB, z1, z2, y1, y2);

                                // Arc match score from LocARNA
                                score_t arc_score = score_t(locarna_scoring_->arcmatch(arc_match));

                                // z-column scoring (sequence layer only)
                                // Structure layer score is in arcmatch_score, not here
                                score_t z_col_score(0);
                                ColumnType col_U_z;
                                score_t u_s_z(0);

                                if (c1_z == 1 && c2_z == 1) {
                                    u_s_z = score_t(locarna_scoring_->basematch(z1, z2));
                                    col_U_z = ColumnType::MATCH;
                                } else if (c1_z == 1 && c2_z == 0) {
                                    u_s_z = score_t(locarna_scoring_->gapA(z1));
                                    col_U_z = ColumnType::DEL_A;
                                } else if (c1_z == 0 && c2_z == 1) {
                                    u_s_z = score_t(locarna_scoring_->gapB(z2));
                                    col_U_z = ColumnType::INS_A;
                                } else {
                                    col_U_z = ColumnType::MATCH;
                                }

                                // Structure layer always has content in Case 2 (MATCH)
                                // but score is 0 here - structure contribution is in arc_score
                                ColumnType col_V_z = ColumnType::MATCH;

                                score_t w_s_z = shift_scoring_.shift_penalty(col_U_z, col_V_z);
                                z_col_score = u_s_z + w_s_z;

                                // y-column scoring (sequence layer only)
                                score_t y_col_score(0);
                                ColumnType col_U_y;
                                score_t u_s_y(0);

                                if (c1_y == 1 && c2_y == 1) {
                                    u_s_y = score_t(locarna_scoring_->basematch(y1, y2));
                                    col_U_y = ColumnType::MATCH;
                                } else if (c1_y == 1 && c2_y == 0) {
                                    u_s_y = score_t(locarna_scoring_->gapA(y1));
                                    col_U_y = ColumnType::DEL_A;
                                } else if (c1_y == 0 && c2_y == 1) {
                                    u_s_y = score_t(locarna_scoring_->gapB(y2));
                                    col_U_y = ColumnType::INS_A;
                                } else {
                                    col_U_y = ColumnType::MATCH;
                                }

                                // Structure layer always has content in Case 2 (MATCH)
                                ColumnType col_V_y = ColumnType::MATCH;

                                score_t w_s_y = shift_scoring_.shift_penalty(col_U_y, col_V_y);
                                y_col_score = u_s_y + w_s_y;

                                // Total score
                                score_t total_score = m_before + d_score + arc_score +
                                                      z_col_score + y_col_score;

                                if (total_score > best_score) {
                                    best_score = total_score;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return best_score;
}

} // namespace RNAShiftAlign
