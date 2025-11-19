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
    // M matrix: 4D matrix for sequence alignment M(y1, y2, y3, y4)
    M_ = std::make_unique<LocARNA::ShiftMatrixM<score_t>>(
        lenA + 1,
        lenB + 1,
        max_shifts_
    );

    // Initialize M matrix with -infinity (maximization problem)
    M_->fill(LocARNA::infty_score_t::neg_infty);
    // Base case: M(0, 0, 0, 0) = 0 (empty alignment)
    M_->set(0, 0, 0, 0, score_t(0));

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
    // Two-phase Sankoff algorithm:
    // 1. Fill D-matrix (processes arcs in descending order of left endpoints)
    // 2. Fill top-level M-matrix using pre-computed D-values

    fill_D();
    fill_M_with_structure();

    // Optimal score is at M(len_A, len_B, len_A, len_B)
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

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

void
ShiftAligner::fill_M_local(size_type left_A, size_type right_A,
                            size_type left_B, size_type right_B) {
    // Fill M-matrix for region [left_A..right_A] x [left_B..right_B]
    // Used when computing D-values for arc matches
    // Optimizes over BOTH unpaired columns AND paired cases (using D-values for inner arcs)

    // Define all 15 valid column types (c1, c2, c3, c4)
    // where c_i ∈ {0, 1} and not all zeros
    // 1 = consume position (•), 0 = gap (-)
    struct Column {
        int c1, c2, c3, c4;  // Column pattern
    };

    // All 15 valid column types (excluding (0,0,0,0))
    // Ordered to prefer no-shift columns (c1==c3 && c2==c4) first
    // This ensures deterministic tie-breaking in favor of U==V
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

    // Fill M matrix for local region
    // Respects dependencies: fill in order of increasing positions

    for (size_type y1 = left_A; y1 <= right_A; ++y1) {
        for (size_type y2 = left_B; y2 <= right_B; ++y2) {
            for (size_type y3 = left_A; y3 <= right_A; ++y3) {
                // Check δ_max constraint
                if (std::abs(static_cast<int>(y1) - static_cast<int>(y3)) > static_cast<int>(max_shifts_))
                    continue;

                for (size_type y4 = left_B; y4 <= right_B; ++y4) {
                    // Check δ_max constraint
                    if (std::abs(static_cast<int>(y2) - static_cast<int>(y4)) > static_cast<int>(max_shifts_))
                        continue;

                    // Skip base case
                    if (y1 == left_A && y2 == left_B && y3 == left_A && y4 == left_B)
                        continue;

                    // Try all valid column types
                    score_t best_score = score_t::neg_infty;

                    for (const auto& col : valid_columns) {
                        // Compute predecessor position x = y - c
                        int x1 = y1 - col.c1;
                        int x2 = y2 - col.c2;
                        int x3 = y3 - col.c3;
                        int x4 = y4 - col.c4;

                        // Check if predecessor is within local region
                        if (x1 < static_cast<int>(left_A) || x2 < static_cast<int>(left_B) ||
                            x3 < static_cast<int>(left_A) || x4 < static_cast<int>(left_B))
                            continue;

                        // Check δ_max constraint at predecessor
                        if (std::abs(x1 - x3) > static_cast<int>(max_shifts_))
                            continue;
                        if (std::abs(x2 - x4) > static_cast<int>(max_shifts_))
                            continue;

                        // Get score from predecessor
                        score_t pred_score = M_->get(x1, x2, x3, x4);

                        // Compute column score s(y, c) = u_s + v_s + w_s
                        score_t column_score(0);

                        // u_s: sequence alignment score for U layer (positions y1, y2)
                        score_t u_s(0);
                        ColumnType col_U;
                        if (col.c1 == 1 && col.c2 == 1) {
                            // Both sequences advance: match or mismatch
                            u_s = score_t(locarna_scoring_->basematch(y1, y2));
                            col_U = ColumnType::MATCH;
                        } else if (col.c1 == 1 && col.c2 == 0) {
                            // Gap in sequence B (deletion from A)
                            u_s = score_t(locarna_scoring_->gapA(y1));
                            col_U = ColumnType::DEL_A;
                        } else if (col.c1 == 0 && col.c2 == 1) {
                            // Gap in sequence A (insertion to B)
                            u_s = score_t(locarna_scoring_->gapB(y2));
                            col_U = ColumnType::INS_A;
                        } else {
                            // Both gaps in U: score 0, column type doesn't matter for shift scoring
                            col_U = ColumnType::MATCH;  // Arbitrary choice for (0,0)
                        }

                        // v_s: structure score for V layer (positions y3, y4)
                        score_t v_s(0);
                        ColumnType col_V;
                        if (col.c3 == 1 && col.c4 == 1) {
                            // Both sequences advance: match or mismatch
                            v_s = score_t(locarna_scoring_->basematch(y3, y4));
                            col_V = ColumnType::MATCH;
                        } else if (col.c3 == 1 && col.c4 == 0) {
                            // Gap in sequence B
                            v_s = score_t(locarna_scoring_->gapA(y3));
                            col_V = ColumnType::DEL_A;
                        } else if (col.c3 == 0 && col.c4 == 1) {
                            // Gap in sequence A
                            v_s = score_t(locarna_scoring_->gapB(y4));
                            col_V = ColumnType::INS_A;
                        } else {
                            // Both gaps in V: score 0
                            col_V = ColumnType::MATCH;  // Arbitrary choice for (0,0)
                        }

                        // w_s: shift penalty based on gap pattern mismatch
                        // Compare column types c_U vs c_V
                        score_t w_s = shift_scoring_.shift_penalty(col_U, col_V);

                        column_score = u_s + v_s + w_s;

                        // Compute total score for this choice
                        score_t total_score = pred_score + column_score;

                        // Update best score
                        if (total_score > best_score) {
                            best_score = total_score;
                        }
                    }

                    // ========== CASE 2: Paired positions ==========
                    // Paired columns: structure layer (y3, y4) must match arc right ends
                    // Valid gap patterns: (1,1,1,1), (1,0,1,1), (0,1,1,1), (0,0,1,1)
                    // Structure layer cannot have gaps when base pairs are matched

                    // Get all arc matches with right ends at current structure positions (y3, y4)
                    const auto& arcmatches_at_right = arc_matches_->common_right_end_list(y3, y4);

                    for (const auto& am_idx : arcmatches_at_right) {
                        const auto& arc_match = arc_matches_->arcmatch(am_idx);
                        const auto& arcA = arc_match.arcA();
                        const auto& arcB = arc_match.arcB();

                        size_type z3 = arcA.left();   // Left arc endpoint in structure A
                        size_type z4 = arcB.left();   // Left arc endpoint in structure B

                        // Arc must be inside current region
                        if (z3 <= left_A || z4 <= left_B)
                            continue;

                        // Iterate over all possible sequence positions z1, z2 at arc left ends
                        // z1, z2 can shift from z3, z4 within δ_max constraint
                        for (int shift_z1 = -static_cast<int>(max_shifts_);
                             shift_z1 <= static_cast<int>(max_shifts_); ++shift_z1) {
                            for (int shift_z2 = -static_cast<int>(max_shifts_);
                                 shift_z2 <= static_cast<int>(max_shifts_); ++shift_z2) {

                                int z1 = static_cast<int>(z3) + shift_z1;
                                int z2 = static_cast<int>(z4) + shift_z2;

                                // Check bounds for z positions
                                if (z1 < static_cast<int>(left_A) || z2 < static_cast<int>(left_B))
                                    continue;

                                // Iterate over valid gap patterns for z-column (arc left ends)
                                // Gap patterns: (c1_z, c2_z, 1, 1) where c1_z, c2_z ∈ {0, 1}
                                for (int c1_z = 0; c1_z <= 1; ++c1_z) {
                                    for (int c2_z = 0; c2_z <= 1; ++c2_z) {

                                        // Iterate over valid gap patterns for y-column (arc right ends)
                                        // Gap patterns: (c1_y, c2_y, 1, 1) where c1_y, c2_y ∈ {0, 1}
                                        for (int c1_y = 0; c1_y <= 1; ++c1_y) {
                                            for (int c2_y = 0; c2_y <= 1; ++c2_y) {

                                                // M score before arc left endpoints: M(z1-c1_z, z2-c2_z, z3-1, z4-1)
                                                int x1 = z1 - c1_z;
                                                int x2 = z2 - c2_z;
                                                int x3 = static_cast<int>(z3) - 1;
                                                int x4 = static_cast<int>(z4) - 1;

                                                // Check bounds for x positions
                                                if (x1 < static_cast<int>(left_A) || x2 < static_cast<int>(left_B) ||
                                                    x3 < static_cast<int>(left_A) || x4 < static_cast<int>(left_B))
                                                    continue;

                                                // Check δ_max constraint at x positions
                                                if (std::abs(x1 - x3) > static_cast<int>(max_shifts_))
                                                    continue;
                                                if (std::abs(x2 - x4) > static_cast<int>(max_shifts_))
                                                    continue;

                                                // Get M score before the arc
                                                score_t m_before = M_->get(x1, x2, x3, x4);

                                                // D score: optimal alignment inside matched arcs
                                                // D indexed by sequence positions at arc endpoints (z1, z2, y1, y2)
                                                score_t d_score = D_->get(arcA, arcB, z1, z2, y1, y2);

                                                // Arc match score from LocARNA (includes structure scoring)
                                                score_t arc_score = score_t(locarna_scoring_->arcmatch(arc_match));

                                                // Score for the two paired columns (z and y)
                                                score_t z_col_score(0);
                                                score_t y_col_score(0);

                                                // z-column scoring (sequence layer)
                                                ColumnType col_U_z, col_V_z;
                                                score_t u_s_z(0), v_s_z(0);

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
                                                    // (0,0) - gap-on-gap scores 0
                                                    col_U_z = ColumnType::MATCH;  // Arbitrary for (0,0)
                                                }

                                                // Structure layer always has content in Case 2
                                                v_s_z = score_t(locarna_scoring_->basematch(z3, z4));
                                                col_V_z = ColumnType::MATCH;

                                                score_t w_s_z = shift_scoring_.shift_penalty(col_U_z, col_V_z);
                                                z_col_score = u_s_z + v_s_z + w_s_z;

                                                // y-column scoring (sequence layer)
                                                ColumnType col_U_y, col_V_y;
                                                score_t u_s_y(0), v_s_y(0);

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
                                                    // (0,0) - gap-on-gap scores 0
                                                    col_U_y = ColumnType::MATCH;  // Arbitrary for (0,0)
                                                }

                                                // Structure layer always has content in Case 2
                                                v_s_y = score_t(locarna_scoring_->basematch(y3, y4));
                                                col_V_y = ColumnType::MATCH;

                                                score_t w_s_y = shift_scoring_.shift_penalty(col_U_y, col_V_y);
                                                y_col_score = u_s_y + v_s_y + w_s_y;

                                                score_t paired_col_score = z_col_score + y_col_score;

                                                // Total score: M(before) + D(inside) + arc_score + column_scores
                                                score_t total_score = m_before + d_score + arc_score + paired_col_score;

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

                    // Store the best score
                    M_->set(y1, y2, y3, y4, best_score);
                }
            }
        }
    }
}

void
ShiftAligner::fill_D() {
    // Implement D matrix computation following LocARNA pattern
    // Process arc matches in descending order of left endpoints
    // For each arc: fill local M, extract D values

    using Arc = LocARNA::BasePairs__Arc;

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Traverse left ends in descending order (inner arcs before outer arcs)
    for (size_type al = lenA; al > 0; --al) {
        for (size_type bl = lenB; bl > 0; --bl) {

            // Get all arc matches with left ends at (al, bl)
            const auto& arcmatches_at_left = arc_matches_->common_left_end_list(al, bl);
            if (arcmatches_at_left.empty())
                continue;

            // Find maximal right ends for arcs starting at (al, bl)
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

            // Fill LOCAL M-matrix for region inside the arc
            // This reuses M_ memory and can use D-values for inner arcs
            fill_M_local(al, max_ar, bl, max_br);

            // Extract D values for all arc matches with left ends at (al, bl)
            for (const auto& am_idx : arcmatches_at_left) {
                const auto& arc_match = arc_matches_->arcmatch(am_idx);
                const Arc& arcA = arc_match.arcA();
                const Arc& arcB = arc_match.arcB();

                size_type right_A = arcA.right();
                size_type right_B = arcB.right();

                // Create offset matrix for this arc pair
                D_->create_offsetmatrix(arcA.idx(), arcB.idx());

                // Extract D values: M at positions just inside arc right ends
                // Structure positions are fixed at (right_A-1, right_B-1)
                // Sequence positions (y1, y2) can vary within δ_max
                size_type y3_fixed = right_A - 1;
                size_type y4_fixed = right_B - 1;

                for (int shift_y1 = -static_cast<int>(max_shifts_);
                     shift_y1 <= static_cast<int>(max_shifts_); ++shift_y1) {
                    for (int shift_y2 = -static_cast<int>(max_shifts_);
                         shift_y2 <= static_cast<int>(max_shifts_); ++shift_y2) {

                        // Sequence positions can shift from structure positions
                        int y1 = static_cast<int>(y3_fixed) + shift_y1;
                        int y2 = static_cast<int>(y4_fixed) + shift_y2;

                        // Bounds checking
                        if (y1 <= static_cast<int>(al) || y2 <= static_cast<int>(bl))
                            continue;
                        if (y1 >= static_cast<int>(right_A) || y2 >= static_cast<int>(right_B))
                            continue;

                        // Get M value from local M-matrix just filled
                        // M(y1, y2, y3_fixed, y4_fixed) = best alignment up to this point
                        score_t m_score = M_->get(y1, y2, y3_fixed, y4_fixed);

                        // Store in D matrix indexed by sequence positions at arc endpoints
                        // D(arcA, arcB, al, bl, y1, y2)
                        D_->set(arcA, arcB, al, bl, y1, y2, m_score);
                    }
                }
            }
        }
    }
}

void
ShiftAligner::fill_M_with_structure() {
    // Fill top-level M-matrix for entire sequence range
    // Same recursion as fill_M_local, but for range [0..lenA] x [0..lenB]

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Simply call fill_M_local with full sequence range
    fill_M_local(0, lenA, 0, lenB);
}

} // namespace RNAShiftAlign
