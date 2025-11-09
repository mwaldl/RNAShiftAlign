#include "shift_aligner.hh"

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

    // Step 5: Allocate shift-aware DP matrix M(y1, y2, y3, y4)
    M_ = std::make_unique<LocARNA::ShiftMatrixM<score_t>>(
        lenA + 1,
        lenB + 1,
        max_shifts_
    );

    // Initialize M matrix with -infinity (maximization problem)
    M_->fill(LocARNA::infty_score_t::neg_infty);
    // Base case: M(0, 0, 0, 0) = 0 (empty alignment)
    M_->set(0, 0, 0, 0, score_t(0));
}

ShiftAligner::score_t
ShiftAligner::align() {
    // Fill M matrix for unpaired-only alignment
    fill_M_unpaired();

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

    const std::vector<Column> valid_columns = {
        {1,1,1,1}, {1,1,1,0}, {1,1,0,1}, {1,1,0,0},
        {1,0,1,1}, {1,0,1,0}, {1,0,0,1}, {1,0,0,0},
        {0,1,1,1}, {0,1,1,0}, {0,1,0,1}, {0,1,0,0},
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

void
ShiftAligner::fill_M_unpaired() {
    // Implement Case 1 of Equation 17 (unpaired positions only)
    // M(x, y) = max over c∈C { M(x, y-c) + s(y, c) }

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Define all 15 valid column types (c1, c2, c3, c4)
    // where c_i ∈ {0, 1} and not all zeros
    // 1 = consume position (•), 0 = gap (-)
    struct Column {
        int c1, c2, c3, c4;  // Column pattern
    };

    // All 15 valid column types (excluding (0,0,0,0))
    const std::vector<Column> valid_columns = {
        {1,1,1,1}, {1,1,1,0}, {1,1,0,1}, {1,1,0,0},
        {1,0,1,1}, {1,0,1,0}, {1,0,0,1}, {1,0,0,0},
        {0,1,1,1}, {0,1,1,0}, {0,1,0,1}, {0,1,0,0},
        {0,0,1,1}, {0,0,1,0}, {0,0,0,1}
    };

    // Fill M matrix in order of increasing positions
    // y = (y1, y2, y3, y4) ranges from (0,0,0,0) to (lenA, lenB, lenA, lenB)
    // We need to fill in an order that respects dependencies

    for (size_type y1 = 0; y1 <= lenA; ++y1) {
        for (size_type y2 = 0; y2 <= lenB; ++y2) {
            for (size_type y3 = 0; y3 <= lenA; ++y3) {
                // Check δ_max constraint: |y1 - y3| ≤ δ_max
                if (std::abs(static_cast<int>(y1) - static_cast<int>(y3)) > static_cast<int>(max_shifts_))
                    continue;

                for (size_type y4 = 0; y4 <= lenB; ++y4) {
                    // Check δ_max constraint: |y2 - y4| ≤ δ_max
                    if (std::abs(static_cast<int>(y2) - static_cast<int>(y4)) > static_cast<int>(max_shifts_))
                        continue;

                    // Skip base case (already initialized)
                    if (y1 == 0 && y2 == 0 && y3 == 0 && y4 == 0)
                        continue;

                    // Try all valid column types
                    score_t best_score = score_t::neg_infty;

                    for (const auto& col : valid_columns) {
                        // Compute predecessor position x = y - c
                        int x1 = y1 - col.c1;
                        int x2 = y2 - col.c2;
                        int x3 = y3 - col.c3;
                        int x4 = y4 - col.c4;

                        // Check if predecessor is valid
                        if (x1 < 0 || x2 < 0 || x3 < 0 || x4 < 0)
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

                    // Store the best score
                    M_->set(y1, y2, y3, y4, best_score);
                }
            }
        }
    }
}

} // namespace RNAShiftAlign
