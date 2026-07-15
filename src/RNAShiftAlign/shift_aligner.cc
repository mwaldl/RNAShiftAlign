/**
 * \file shift_aligner.cc
 *
 * \brief Implementation of the Sankoff-style RNA bi-aligner.
 *
 * Layout of this file:
 *  - constructor: preprocessing (fold both RNAs, enumerate arc matches, build
 *    the LocARNA scoring scheme, allocate D)
 *  - align(): the two-phase forward pass
 *  - traceback(), trace_M(), trace_M_paired(), trace_D(), append_column():
 *    reconstruction of the U/V layers and the consensus structure
 *  - align_D(), align_in_arcmatch(), fill_D_entries(): filling the DP matrices
 *  - compute_unpaired_score(), compute_paired_score(): the two recursion cases
 *
 * Coordinates used throughout: (y1, y2) are sequence-layer positions in A and B,
 * (y3, y4) are structure-layer positions, and (al, ar, bl, br) are the structure
 * boundaries of the arc region being filled. See shift_aligner.hh for the column
 * model and the δ_max shift bound.
 * 
 * Implementation builds on Locarna library by Sebastian Will.
 *
 * Copyright (C) Maria Waldl <code@waldl.org>
 */

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
                           const ScoringParams &scoring_params,
                           int verbose)
    : shift_scoring_(score_t(scoring_params.delta)),
      max_shifts_(scoring_params.max_shifts),
      verbose_(verbose),
      alignment_score_(score_t(0))
{
    if (verbose_ >= 1) {
        std::cout << "=== ShiftAligner Constructor ===" << std::endl;
    }

    // Create sequences
    seqA_ = std::make_unique<LocARNA::Sequence>("seqA", seqA);
    seqB_ = std::make_unique<LocARNA::Sequence>("seqB", seqB);

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    if (verbose_ >= 2) {
        std::cout << "Sequence lengths: A=" << lenA << ", B=" << lenB << std::endl;
    }

    // Compute base pair probabilities via ViennaRNA partition function
    // Sequence → MultipleAlignment → RnaEnsemble → RnaData
    // Note: using the multiple alignment repesentation is currently not needed, 
    //       but enables future extension to multiple sequence alignment.
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

    // Enumerate arc matches between the two sequences.
    // ArcMatches requires an AnchorConstraints and a TraceController; both are
    // constructed unrestricted here, since neither anchoring nor banding is
    // exposed by this tool.
    LocARNA::AnchorConstraints seq_constraints(lenA, "", lenB, "", true);

    auto trace_controller_ = std::make_unique<LocARNA::TraceController>(
        *seqA_,
        *seqB_,
        nullptr,  // no reference alignment
        -1,       // max_diff = -1: no banding
        false     // no relaxation
    );

    // TODO: implement use of tracecontroler
    //TraceController trace_controller(seqA, seqB, multiple_ref_alignment.get(),
    //                                 clp.max_diff, clp.max_diff_relax);
    //restrict_trace_by_probabilities(clp, rna_dataA.get(), rna_dataB.get(),
    //                                ribosum.get(), ribofit.get(),
    //                                &trace_controller,
    //                                pf_score_t());


    // Arc-match filtering is controlled by min_prob and the max_diff_* limits
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
        *trace_controller_,
        seq_constraints
    );

    if (verbose_ >= 3) {
        std::cout << "\n=== Base Pairs and Arc Matches (verbose level 3) ===" << std::endl;

        // Output base pairs for sequence A
        const auto& bps_A = arc_matches_->get_base_pairsA();
        std::cout << "\nBase pairs in sequence A (probability > " << prep_params.min_prob << "):" << std::endl;
        std::cout << "Total base pairs: " << bps_A.num_bps() << std::endl;
        for (size_type i = 0; i < bps_A.num_bps(); ++i) {
            const auto& arc = bps_A.arc(i);
            double prob = rna_dataA_->arc_prob(arc.left(), arc.right());
            std::cout << "  BP " << i << ": (" << arc.left() << ", " << arc.right()
                      << ") prob=" << prob << std::endl;
        }

        // Output base pairs for sequence B
        const auto& bps_B = arc_matches_->get_base_pairsB();
        std::cout << "\nBase pairs in sequence B (probability > " << prep_params.min_prob << "):" << std::endl;
        std::cout << "Total base pairs: " << bps_B.num_bps() << std::endl;
        for (size_type i = 0; i < bps_B.num_bps(); ++i) {
            const auto& arc = bps_B.arc(i);
            double prob = rna_dataB_->arc_prob(arc.left(), arc.right());
            std::cout << "  BP " << i << ": (" << arc.left() << ", " << arc.right()
                      << ") prob=" << prob << std::endl;
        }

        // Output arc matches
        std::cout << "\nArc matches between sequences:" << std::endl;
        std::cout << "Total arc matches: " << arc_matches_->num_arc_matches() << std::endl;
        for (size_type i = 0; i < arc_matches_->num_arc_matches(); ++i) {
            const LocARNA::ArcMatch& am = arc_matches_->arcmatch(i);
            const LocARNA::BasePairs__Arc& arcA = am.arcA();
            const LocARNA::BasePairs__Arc& arcB = am.arcB();
            std::cout << "  Match " << i << ": A(" << arcA.left() << ", " << arcA.right()
                      << ") <-> B(" << arcB.left() << ", " << arcB.right() << ")" << std::endl;
        }
        std::cout << std::endl;
    }

    // Initialize scoring scheme using named argument pattern
    // Compute expected probabilities
    double exp_probA = rna_dataA_->arc_cutoff_prob();
    double exp_probB = rna_dataB_->arc_cutoff_prob();

    // Create scoring parameters using named argument pattern.
    // Must be stored as a member: LocARNA::Scoring keeps only a pointer to
    // this object, so it must outlive locarna_scoring_.
    locarna_scoring_params_ = std::make_unique<LocARNA::ScoringParams>(
        LocARNA::ScoringParams::match(scoring_params.match),
        LocARNA::ScoringParams::mismatch(scoring_params.mismatch),
        LocARNA::ScoringParams::indel(scoring_params.indel),
        //LocARNA::ScoringParams::indel_opening(scoring_params.indel_opening),
        LocARNA::ScoringParams::struct_weight(scoring_params.struct_weight),
        //LocARNA::ScoringParams::tau_factor(scoring_params.tau_factor),
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
        *locarna_scoring_params_
    );

    // Allocate shift-aware DP matrices
    // M matrix is allocated fresh in each align_in_arcmatch() call
    // since entries are independent between different starting columns

    // D matrix: 6D matrix for structure alignment D(arc_a, arc_b, x1, x2, y1, y2)
    // Dimensions: number of base pairs in each RNA
    size_type num_bps_A = arc_matches_->get_base_pairsA().num_bps();
    size_type num_bps_B = arc_matches_->get_base_pairsB().num_bps();
    D_ = std::make_unique<ShiftMatrixD<score_t>>(
        num_bps_A,
        num_bps_B,
        max_shifts_
    );
    // Note: Individual 4D offset matrices are created on-demand in fill_D()

    if (verbose_ >= 2)
        std::cout << "Preprocessing complete: "
                  << arc_matches_->num_arc_matches() << " arc matches" << std::endl;
}

ShiftAligner::score_t
ShiftAligner::align() {
    // Two-phase Sankoff algorithm following LocARNA pattern:
    // 1. Fill D-matrix (processes arcs in descending order of left endpoints)
    // 2. Align top level with imagined surrounding arc (0, lenA+1, 0, lenB+1)

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Phase 1: Fill D matrix for all arc matches
    if (verbose_ >= 2)
        std::cout << "Phase 1: filling D matrix" << std::endl;

    align_D();

    if (verbose_ >= 3) {
        std::cout << "\n=== D matrix after align_D ===" << std::endl;
        D_->debug_print(std::cout, [](const score_t& s) { return s == score_t::neg_infty; });
    }

    // Phase 2: Top-level alignment with imagined arc around entire sequences
    // Structure boundaries: (0, lenA+1, 0, lenB+1) - virtual arc endpoints
    // Sequence interior: (1, 1) to (lenA, lenB) - actual sequence positions
    if (verbose_ >= 2)
        std::cout << "Phase 2: top-level alignment" << std::endl;

    align_in_arcmatch(0, lenA + 1, 0, lenB + 1, // pseudo arch match enclosing the alignment
                      1, 1,           // x1, x2: first interior sequence positions
                      lenA, lenB);    // y1, y2: sequence end

    // Optimal score is at M(lenA, lenB, lenA, lenB)
    // Both sequence and structure layers must consume all positions
    alignment_score_ = M_->get(lenA, lenB, lenA, lenB);
    return alignment_score_;
}

// ============================================================================
// TRACEBACK IMPLEMENTATION (LocARNA-style recursive structure)
// ============================================================================

void
ShiftAligner::traceback() {
    // Clear previous alignments
    alignment_U_seqA_.clear();
    alignment_U_seqB_.clear();
    alignment_V_seqA_.clear();
    alignment_V_seqB_.clear();
    consensus_structure_.clear();
    dot_bracket_A_.clear();
    dot_bracket_B_.clear();
    arc_stack_.clear();

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Precompute per-position bracket characters (1-based) for prob >= 0.5
    bracket_A_.assign(lenA + 1, '.');
    bracket_B_.assign(lenB + 1, '.');
    {
        const LocARNA::BasePairs& bpsA = arc_matches_->get_base_pairsA();
        for (size_type i = 0; i < bpsA.num_bps(); ++i) {
            const LocARNA::BasePairs__Arc& arc = bpsA.arc(i);
            if (rna_dataA_->arc_prob(arc.left(), arc.right()) >= 0.5) {
                bracket_A_[arc.left()]  = '(';
                bracket_A_[arc.right()] = ')';
            }
        }
    }
    {
        const LocARNA::BasePairs& bpsB = arc_matches_->get_base_pairsB();
        for (size_type i = 0; i < bpsB.num_bps(); ++i) {
            const LocARNA::BasePairs__Arc& arc = bpsB.arc(i);
            if (rna_dataB_->arc_prob(arc.left(), arc.right()) >= 0.5) {
                bracket_B_[arc.left()]  = '(';
                bracket_B_[arc.right()] = ')';
            }
        }
    }


    if (verbose_ >= 2) {
        std::cout << "\n=== Starting Traceback ===" << std::endl;
        std::cout << "Endpoint: M(" << lenA << "," << lenB << "," << lenA << "," << lenB << ")" << std::endl;
        std::cout << "Score: " << M_->get(lenA, lenB, lenA, lenB) << std::endl;
    }

    // Start traceback from endpoint: M(lenA, lenB, lenA, lenB)
    // Structure boundaries: (0, lenA+1, 0, lenB+1) - imagined surrounding arc
    // Sequence boundaries: (0, 0) to (lenA, lenB)
    trace_M(0, lenA + 1,    // al, ar: structure boundaries for A
            0, lenB + 1,    // bl, br: structure boundaries for B
            1, 1,           // x1, x2: sequence starting positions
            lenA, lenB,     // y1, y2: current sequence positions
            lenA, lenB);    // y3, y4: current structure positions

    if (verbose_ >= 2) {
        std::cout << "Traceback completed. Alignment length: " << alignment_U_seqA_.length() << std::endl;
    }
}

void
ShiftAligner::trace_M(size_type al, size_type ar,
                      size_type bl, size_type br,
                      size_type x1, size_type x2,
                      size_type y1, size_type y2,
                      size_type y3, size_type y4) {
    // Base case: reached starting position.
    // M is initialized at (x1-1, x2-1, al, bl) = 0 in init_M.
    if (y1 == x1 - 1 && y2 == x2 - 1 && y3 == al && y4 == bl) {
        return;
    }

    if (verbose_ >= 3) {
        std::cout << "trace_M: (" << y1 << "," << y2 << "," << y3 << "," << y4 << ")" << std::endl;
    }

    score_t current_score = M_->get(y1, y2, y3, y4);

    // Define column types (same ordering as forward pass)
    struct Column {
        int c1, c2, c3, c4;
    };

    const std::vector<Column> valid_columns = {
        // No-shift columns first
        {1,1,1,1}, {1,0,1,0}, {0,1,0,1},
        // Shift columns
        {1,1,1,0}, {1,1,0,1}, {1,1,0,0},
        {1,0,1,1}, {1,0,0,1}, {1,0,0,0},
        {0,1,1,1}, {0,1,1,0}, {0,1,0,0},
        {0,0,1,1}, {0,0,1,0}, {0,0,0,1}
    };

    // Case 1: Try all unpaired column types
    for (const auto& col : valid_columns) {
        // Compute predecessor
        int px1 = static_cast<int>(y1) - col.c1;
        int px2 = static_cast<int>(y2) - col.c2;
        int px3 = static_cast<int>(y3) - col.c3;
        int px4 = static_cast<int>(y4) - col.c4;

        // Validate predecessor.
        if (px1 < static_cast<int>(x1) - 1 || px2 < static_cast<int>(x2) - 1 ||
            px3 < static_cast<int>(al) || px4 < static_cast<int>(bl))
            continue;



        // Check δ_max constraint
        if (std::abs(px1 - px3) > static_cast<int>(max_shifts_))
            continue;
        if (std::abs(px2 - px4) > static_cast<int>(max_shifts_))
            continue;

        // Get predecessor score
        score_t pred_score = M_->get(px1, px2, px3, px4);

        // Compute column score
        score_t u_s(0);
        if (col.c1 == 1 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->basematch(y1, y2));
        } else if (col.c1 == 1 && col.c2 == 0) {
            u_s = score_t(locarna_scoring_->gapA(y1));
        } else if (col.c1 == 0 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->gapB(y2));
        }
        // else: (0,0) gap-on-gap, u_s = 0

        score_t v_s(0);
        if (col.c3 == 1 && col.c4 == 1) {
            // v_s = score_t(locarna_scoring_->basematch(y3, y4)); // 0 for match/mismatch
        } else if (col.c3 == 1 && col.c4 == 0) {
            v_s = score_t(locarna_scoring_->gapA(y3));
        } else if (col.c3 == 0 && col.c4 == 1) {
            v_s = score_t(locarna_scoring_->gapB(y4));
        }
        // else: (0,0) gap-on-gap, v_s = 0

        // Compute shift penalty using column components directly
        score_t w_s = shift_scoring_.shift_penalty(col.c1, col.c2, col.c3, col.c4);
        score_t column_score = u_s + v_s + w_s;

        // Check if this is the correct predecessor
        if (pred_score + column_score == current_score) {
            // Found it! Recurse on predecessor
            trace_M(al, ar, bl, br, x1, x2, px1, px2, px3, px4);

            // Append column to alignment (building forward)
            append_column(col.c1, col.c2, col.c3, col.c4, y1, y2, y3, y4);
            return;
        }
    }

    if (verbose_ >= 3) {
        std::cout << "trace_M: no matching unpaired case" << std::endl;
    }

    // Case 2: Try paired/arc match cases
    trace_M_paired(al, ar, bl, br, x1, x2, y1, y2, y3, y4, current_score);
}

void
ShiftAligner::trace_M_paired(size_type al, size_type ar,
                             size_type bl, size_type br,
                             size_type x1, size_type x2,
                             size_type y1, size_type y2,
                             size_type y3, size_type y4,
                             score_t current_score) {
    // Get arc matches with right ends at (y3, y4)
    const LocARNA::ArcMatchIdxVec& arcmatch_right_idx =
        arc_matches_->common_right_end_list(y3, y4);

    for (LocARNA::ArcMatch::idx_type am_idx : arcmatch_right_idx) {
        const LocARNA::ArcMatch& am = arc_matches_->arcmatch(am_idx);
        const LocARNA::BasePairs__Arc& arcA = am.arcA();
        const LocARNA::BasePairs__Arc& arcB = am.arcB();

        size_type z3 = arcA.left();
        size_type z4 = arcB.left();

        if (verbose_ >= 3)
            std::cout << "trace_M_paired: trying arcmatch " << am_idx
                      << " " << arcA << " " << arcB << std::endl;

        // Arc must be inside current region
        if (z3 <= al || z4 <= bl)
            continue;

        // Try all valid sequence positions (z1, z2) at arc left ends
        // D matrix already optimized over all gap pattern combinations
        for (int shift_z1 = -static_cast<int>(max_shifts_);
             shift_z1 <= static_cast<int>(max_shifts_); ++shift_z1) {
            for (int shift_z2 = -static_cast<int>(max_shifts_);
                 shift_z2 <= static_cast<int>(max_shifts_); ++shift_z2) {

                int z1 = static_cast<int>(z3) + shift_z1;
                int z2 = static_cast<int>(z4) + shift_z2;

                // Bounds check: predecessor at (z1-1, z2-1) must be valid
                if (z1 < static_cast<int>(x1) || z2 < static_cast<int>(x2))
                    continue;
        

                // Get D matrix entry - contains complete arc score
                // (left column + interior + right column + arc match)
                score_t d_score = D_->get(arcA, arcB, z1, z2, y1, y2);
                if (d_score == score_t::neg_infty)
                    continue;

                // Get predecessor M score at (z1-1, z2-1, z3-1, z4-1)
                score_t m_before = M_->get(z1 - 1, z2 - 1, z3 - 1, z4 - 1);

                // Check if this path matches current score
                score_t total_score = m_before + d_score;

                if (total_score == current_score) {
                    // Found it. Trace the predecessor first (unless it is the
                    // base case), then let trace_D emit the arc columns.
                    if (verbose_ >= 3)
                        std::cout << "trace_M_paired: took arcmatch " << am_idx << std::endl;

                    if (!(z1 - 1 < static_cast<int>(x1) && z2 - 1 < static_cast<int>(x2) &&
                          z3 - 1 == static_cast<int>(al) && z4 - 1 == static_cast<int>(bl))) {
                        trace_M(al, ar, bl, br, x1, x2, z1 - 1, z2 - 1, z3 - 1, z4 - 1);
                    }

                    // trace_D reconstructs gap patterns and adds arc endpoint columns
                    trace_D(arcA, arcB, z1, z2, y1, y2,
                            score_t(locarna_scoring_->arcmatch(am)));
                    return;
                }
            }
        }
    }

    // No valid traceback found
    throw std::runtime_error("Traceback failed: no valid arc match found at (" +
                            std::to_string(y1) + "," + std::to_string(y2) + "," +
                            std::to_string(y3) + "," + std::to_string(y4) + ")");
}

void
ShiftAligner::trace_D(const LocARNA::BasePairs__Arc& arcA,
                      const LocARNA::BasePairs__Arc& arcB,
                      size_type x1, size_type x2,
                      size_type y1, size_type y2,
                      score_t arc_match_score) {
    // Structure positions
    size_type al = arcA.left();
    size_type ar = arcA.right();
    size_type bl = arcB.left();
    size_type br = arcB.right();

    if (verbose_ >= 3) {
        std::cout << "trace_D: arcA(" << al << "," << ar << ") arcB("
                  << bl << "," << br << ") seq(" << x1 << "," << x2
                  << ")->(" << y1 << "," << y2 << ")" << std::endl;
    }

    // Only D is retained from the forward pass, so this arc's local M must be
    // rebuilt before the interior can be traced. The interior's first sequence
    // position depends on the arc's left gap pattern (a consumed base advances
    // it, a gap does not), and that pattern is not stored — so both it and the
    // right pattern are recovered below by finding the combination that
    // reproduces the D value the forward pass recorded.
    //
    // The y extents must match align_D exactly, or the rebuilt M will differ
    // from the one D was derived from.
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();
    size_type max_y1 = std::min((ar - 1) + max_shifts_, lenA);
    size_type max_y2 = std::min((br - 1) + max_shifts_, lenB);

    score_t d_score_stored = D_->get(arcA, arcB, x1, x2, y1, y2);

    if (verbose_ >= 3)
        std::cout << "trace_D: stored D value: " << d_score_stored << std::endl;

    size_type ar_interior = ar - 1;
    size_type br_interior = br - 1;

    int best_c1_left = 1, best_c2_left = 1, best_c1_right = 1, best_c2_right = 1;

    // Try all 4 left gap patterns, starting with most likely (1,1) = both bases present
    // Order: (1,1), (1,0), (0,1), (0,0) - gaps are rare, so check no-gap case first
    for (int c1_left = 1; c1_left >= 0; --c1_left) {
        for (int c2_left = 1; c2_left >= 0; --c2_left) {

            // Compute interior starting positions based on left gap pattern
            size_type x1_interior = x1 + c1_left;  // If c1_left=1, x1 consumed, interior at x1+1
            size_type x2_interior = x2 + c2_left;  // If c2_left=1, x2 consumed, interior at x2+1

            // Recreate M matrix with correct interior starting positions
            align_in_arcmatch(al, ar, bl, br, x1_interior, x2_interior, max_y1, max_y2);

            // Try all 4 right gap patterns, starting with most likely (1,1) = both bases present
            for (int c1_right = 1; c1_right >= 0; --c1_right) {
                for (int c2_right = 1; c2_right >= 0; --c2_right) {

                    // Compute M predecessor based on right gap pattern
                    int m_y1 = (c1_right == 1) ? static_cast<int>(y1) - 1 : static_cast<int>(y1);
                    int m_y2 = (c2_right == 1) ? static_cast<int>(y2) - 1 : static_cast<int>(y2);

                    // Bounds check: M matrix ranges from (x1_interior-1, x2_interior-1) to (max_y1, max_y2)
                    // Lower bound: must be >= base case position
                    if (m_y1 < static_cast<int>(x1_interior) - 1 || m_y2 < static_cast<int>(x2_interior) - 1)
                        continue;
                    // Upper bound: must be <= what was filled
                    if (m_y1 > static_cast<int>(max_y1) || m_y2 > static_cast<int>(max_y2))
                        continue;

                    score_t m_score = M_->get(m_y1, m_y2, ar_interior, br_interior);
                    if (m_score == score_t::neg_infty)
                        continue;

                    // Compute left column score (sequence positions are always x1, x2)
                    score_t left_seq(0);
                    if (c1_left == 1 && c2_left == 1)
                        left_seq = score_t(locarna_scoring_->basematch(x1, x2));
                    else if (c1_left == 1)
                        left_seq = score_t(locarna_scoring_->gapA(x1));
                    else if (c2_left == 1)
                        left_seq = score_t(locarna_scoring_->gapB(x2));

                    score_t left_shift = shift_scoring_.shift_penalty(c1_left, c2_left, 1, 1);

                    // Compute right column score
                    score_t right_seq(0);
                    if (c1_right == 1 && c2_right == 1)
                        right_seq = score_t(locarna_scoring_->basematch(y1, y2));
                    else if (c1_right == 1)
                        right_seq = score_t(locarna_scoring_->gapA(y1));
                    else if (c2_right == 1)
                        right_seq = score_t(locarna_scoring_->gapB(y2));

                    score_t right_shift = shift_scoring_.shift_penalty(c1_right, c2_right, 1, 1);

                    // Total score
                    score_t total = left_seq + left_shift + m_score + right_seq + right_shift + arc_match_score;

                    if (total == d_score_stored) {
                        // Found it!
                        best_c1_left = c1_left;
                        best_c2_left = c2_left;
                        best_c1_right = c1_right;
                        best_c2_right = c2_right;
                        goto found_pattern;
                    }
                }
            }
        }
    }

    throw std::runtime_error("trace_D: Failed to reconstruct gap patterns for arc");

found_pattern:
    if (verbose_ >= 3) {
        std::cout << "trace_D: Gap patterns: left=(" << best_c1_left << "," << best_c2_left
                  << ") right=(" << best_c1_right << "," << best_c2_right << ")" << std::endl;
    }

    // Recompute interior starting positions with best gap pattern
    size_type x1_interior = x1 + best_c1_left;
    size_type x2_interior = x2 + best_c2_left;

    // Add left arc column
    append_column(best_c1_left, best_c2_left, 1, 1, x1, x2, al, bl, true, false);

    // Trace inside the arc to the specific (y1, y2) position
    // Only trace if not at the base case M(0,0,0,0)
    int m_y1 = (best_c1_right == 1) ? static_cast<int>(y1) - 1 : static_cast<int>(y1);
    int m_y2 = (best_c2_right == 1) ? static_cast<int>(y2) - 1 : static_cast<int>(y2);
    size_type m_y3 = ar - 1;  // Interior structure right
    size_type m_y4 = br - 1;

    if (!(m_y1 < static_cast<int>(x1_interior) && m_y2 < static_cast<int>(x2_interior) &&
          m_y3 == al && m_y4 == bl)) {
        // Pass original structure boundaries (al, ar, bl, br) like in forward pass
        trace_M(al, ar, bl, br,
                x1_interior, x2_interior, m_y1, m_y2, m_y3, m_y4);
    }

    // Add right arc column
    append_column(best_c1_right, best_c2_right, 1, 1, y1, y2, ar, br, false, true);
}

void
ShiftAligner::append_column(int c1, int c2, int c3, int c4,
                            size_type y1, size_type y2,
                            size_type y3, size_type y4,
                            bool is_arc_left,
                            bool is_arc_right) {
    // U layer (sequence alignment)
    if (c1 == 1) {
        alignment_U_seqA_ += seqA_->seqentry(0).seq()[y1];
    } else {
        alignment_U_seqA_ += '-';
    }

    if (c2 == 1) {
        alignment_U_seqB_ += seqB_->seqentry(0).seq()[y2];
    } else {
        alignment_U_seqB_ += '-';
    }

    // V layer (structure alignment)
    if (c3 == 1) {
        alignment_V_seqA_ += seqA_->seqentry(0).seq()[y3];
    } else {
        alignment_V_seqA_ += '-';
    }

    if (c4 == 1) {
        alignment_V_seqB_ += seqB_->seqentry(0).seq()[y4];
    } else {
        alignment_V_seqB_ += '-';
    }

    // Consensus structure annotation
    if (is_arc_left) {
        consensus_structure_ += '(';
    } else if (is_arc_right) {
        consensus_structure_ += ')';
    } else {
        consensus_structure_ += '.';
    }

    // Dot-bracket for each sequence, aligned to structure layer
    dot_bracket_A_ += (c3 == 1) ? bracket_A_[y3] : '-';
    dot_bracket_B_ += (c4 == 1) ? bracket_B_[y4] : '-';
}

std::string
ShiftAligner::format_alignment() const {
    std::ostringstream oss;
    const std::string shifts        = get_shift_string();
    const std::string consensus_seq = get_consensus_sequence();

    // Format 4-way alignment output
    oss << "4-Way Bi-Alignment:\n";
    oss << "U_A:  " << alignment_U_seqA_ << "\n";
    oss << "Cseq: " << consensus_seq << "\n";
    oss << "U_B:  " << alignment_U_seqB_ << "\n";
    oss << "Shft: " << shifts << "\n";
    oss << "V_A:  " << alignment_V_seqA_ << "\n";
    oss << "db_A: " << dot_bracket_A_ << "\n";
    oss << "Cons: " << consensus_structure_ << "\n";
    oss << "db_B: " << dot_bracket_B_ << "\n";
    oss << "V_B:  " << alignment_V_seqB_ << "\n";
    oss << "\n";
    oss << "Legend:\n";
    oss << "  U_A/U_B: Sequence alignment layer (U)\n";
    oss << "  Cseq: consensus sequence (nucleotide if match, '.' if mismatch/gap), based on U\n";
    oss << "  V_A/V_B: Structure alignment layer (V)\n";
    oss << "  db_A/db_B: Dot-bracket representation of each RNA ensemble as centroid (base pairs with prob >= 0.5), aligned to V layer\n";
    oss << "  Cons: Consensus structure ('(' and ')' for paired, '.' for unpaired)\n";
    oss << "  Shft: '=' indicates U==V (no shift), 'S' indicates shift\n";

    return oss.str();
}





// ============================================================================
// Forward pass
//
// Arc endpoints have fixed structure positions (al, bl, ar, br), but the
// sequence positions aligned to them (x1, x2, y1, y2) may drift by up to δ_max.
// The forward pass therefore enumerates those endpoint offsets explicitly,
// which is the main difference from LocARNA's Sankoff recursion.
// ============================================================================

void
ShiftAligner::align_D() {
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Iterate over descending left endpoints: an arc's interior is complete in D before any
    // enclosing arc reads it.
    for (size_type al = lenA; al > 0; --al) {
        // TODO: implement usage of tracecontroler
        //size_type max_bl = std::min(lenB, trace_controller_->max_col(al));
        //size_type min_bl = std::max(1, trace_controller_->min_col(al));
        //for (size_type bl = max_bl; bl >= min_bl; --bl) {
        for (size_type bl = lenB; bl > 0; --bl) {


            // Get arc matches with left ends at (al, bl)
            const LocARNA::ArcMatchIdxVec& arcmatches_at_left =
                arc_matches_->common_left_end_list(al, bl);
            if (arcmatches_at_left.empty())
                continue;
            //if (!(params_->trace_controller_->is_valid_match(al, bl)))
            //    continue;

            // The M matrix is filled once, out to the furthest right end of any
            // arc starting here; shorter arcs read their own right ends from it.
            size_type max_ar = al;
            size_type max_br = bl;
            for (LocARNA::ArcMatch::idx_type am_idx : arcmatches_at_left) {
                const LocARNA::ArcMatch& am = arc_matches_->arcmatch(am_idx);
                max_ar = std::max(max_ar, am.arcA().right());
                max_br = std::max(max_br, am.arcB().right());
            }

            // Degenerate: no arc actually spans an interior here -> skip
            if (max_ar == al || max_br == bl)
                continue;

            // Iterate over sequence offsets at left end (x1, x2)
            // x1, x2 are computed as shifts relative to interior structure positions (al+1, bl+1)
            for (int shift_x1 = -static_cast<int>(max_shifts_);
                 shift_x1 <= static_cast<int>(max_shifts_); ++shift_x1) {
                for (int shift_x2 = -static_cast<int>(max_shifts_);
                     shift_x2 <= static_cast<int>(max_shifts_); ++shift_x2) {

                    // Compute x1, x2 as shifts relative to interior structure start (al+1, bl+1)
                    // This ensures |x1 - (al+1)| <= maxshift for M matrix constraint
                    int x1 = static_cast<int>(al + 1) + shift_x1;
                    int x2 = static_cast<int>(bl + 1) + shift_x2;

                    // Bounds check for x1, x2
                    if (x1 < 1 || x2 < 1)
                        continue;
                    if (x1 > static_cast<int>(lenA) || x2 > static_cast<int>(lenB))
                        continue;

                    // Compute maximum y1, y2 (right end sequence positions)
                    // M matrix interior ends at (ar-1, br-1), so max sequence positions
                    // must satisfy shift constraint relative to these interior boundaries
                    size_type max_y1 = std::min((max_ar - 1) + max_shifts_, lenA);
                    size_type max_y2 = std::min((max_br - 1) + max_shifts_, lenB);

                    if (verbose_ >= 2)
                        std::cout << "  align_D: al=" << al << " bl=" << bl
                                  << " x1=" << x1 << " x2=" << x2 << std::endl;

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
                                 size_type max_y1, size_type max_y2) {
    // Fill M for the region enclosed by the arc match (al, ar) x (bl, br), with
    // the sequence layer starting at (x1, x2).
    //
    //   structure range: y3 in [al, ar-1],   y4 in [bl, br-1]
    //   sequence range:  y1 in [x1-1, max_y1], y2 in [x2-1, max_y2]
    //   base case:       M(x1-1, x2-1, al, bl) = 0  (nothing consumed yet)
    //
    // A single ascending (y1, y2, y3, y4) loop fills every cell, boundary faces
    // included. Two properties make that correct:
    //
    //  - Coverage: clamping y3 to [al, ar-1] intersected with [y1-δ, y1+δ] to
    //    get the set of y3 that are both inside the arc and within δ_max of
    //    y1.
    //  - Ordering: a predecessor y-c is <= y in every coordinate, so ascending
    //    order always computes it first.

    if (verbose_ >= 3) {
        std::cout << "    align_in_arcmatch: arc A(" << al << "," << ar << ") B(" << bl << "," << br << ")"
                  << " x=(" << x1 << "," << x2 << ") max_y=(" << max_y1 << "," << max_y2 << ")" << std::endl;
    }

    // Create fresh M matrix for this arc region. TODO: only allocate slices.
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    M_ = std::make_unique<ShiftMatrixM<score_t>>(lenA + 1, lenB + 1, max_shifts_);
    M_->fill(LocARNA::infty_score_t::neg_infty);

    // Base case: no interior consumed yet
    const size_type x3 = al + 1;   // first interior structure position A
    const size_type x4 = bl + 1;   // first interior structure position B

    M_->set(x1 - 1, x2 - 1, al, bl, score_t(0));

    // Fill M: process all cells in ascending (y1, y2, y3, y4) order.
    const int delta = static_cast<int>(max_shifts_);

    for (size_type y1 = x1 - 1; y1 <= max_y1; ++y1) {
        for (size_type y2 = x2 - 1; y2 <= max_y2; ++y2) {

            // y3 range: [al, ar-1] clamped to shift constraint [y1-δ, y1+δ]
            int y3_lo = std::max(static_cast<int>(al),
                                 static_cast<int>(y1) - delta);
            int y3_hi = std::min(static_cast<int>(ar) - 1,
                                 static_cast<int>(y1) + delta);

            for (int iy3 = y3_lo; iy3 <= y3_hi; ++iy3) {
                size_type y3 = static_cast<size_type>(iy3);

                // y4 range: [bl, br-1] clamped to shift constraint [y2-δ, y2+δ]
                int y4_lo = std::max(static_cast<int>(bl),
                                     static_cast<int>(y2) - delta);
                int y4_hi = std::min(static_cast<int>(br) - 1,
                                     static_cast<int>(y2) + delta);

                for (int iy4 = y4_lo; iy4 <= y4_hi; ++iy4) {
                    size_type y4 = static_cast<size_type>(iy4);

                    // Skip base case — already set to 0
                    if (y1 == x1 - 1 && y2 == x2 - 1 && y3 == al && y4 == bl)
                        continue;

                    // Evaluate both unpaired (Case 1) and paired (Case 2)
                    score_t score = std::max(
                        compute_unpaired_score(x1, x2, x3, x4, y1, y2, y3, y4),
                        compute_paired_score(x1, x2, x3, x4, y1, y2, y3, y4));

                    M_->set(y1, y2, y3, y4, score);
                }
            }
        }
    }
}

void
ShiftAligner::fill_D_entries(size_type al, size_type bl,
                              size_type x1, size_type x2) {
    // Read D entries out of the M matrix that align_D just filled for (al, bl).
    // For every arc match starting here, enumerate the left gap pattern, the
    // right endpoint offsets, and the right gap pattern, and record
    //
    //   D = left column + M(interior right boundary) + right column + arc match
    //
    // keeping the best score per offset key. Several combinations can map to the
    // same key, and the same key is revisited across align_D's (x1, x2) sweep,
    // hence the max-with-existing rather than a plain assignment.

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    const LocARNA::ArcMatchIdxVec& arcmatches_at_left =
        arc_matches_->common_left_end_list(al, bl);

    for (LocARNA::ArcMatch::idx_type am_idx : arcmatches_at_left) {
        const LocARNA::ArcMatch& am = arc_matches_->arcmatch(am_idx);
        const LocARNA::BasePairs__Arc& arcA = am.arcA();
        const LocARNA::BasePairs__Arc& arcB = am.arcB();

        size_type ar = arcA.right();
        size_type br = arcB.right();

        // Create offset matrix for this arc pair if not already allocated
        if (!D_->is_offsetmatrix_allocated(arcA.idx(), arcB.idx())) {
            D_->create_offsetmatrix(arcA.idx(), arcB.idx());
        }

        // Arc match score (base pair structure probabilities) - computed once
        score_t arc_score = score_t(locarna_scoring_->arcmatch(am));

        if (verbose_ >= 3) {
            std::cout << "fill_D_entries: am_idx=" << am_idx
                      << " arcA(" << arcA.left() << "," << ar << ")"
                      << " arcB(" << arcB.left() << "," << br << ")"
                      << " arc_score=" << arc_score << std::endl;
        }

        // Interior structure boundaries
        size_type ar_interior = ar - 1;
        size_type br_interior = br - 1;

        // Iterate over left column gap patterns (4 patterns)
        // Pattern determines which sequence positions are consumed at left endpoint
        for (int c1_left = 0; c1_left <= 1; ++c1_left) {
            for (int c2_left = 0; c2_left <= 1; ++c2_left) {

                int left_seq1 = (c1_left == 1) ? static_cast<int>(x1) - 1 : static_cast<int>(x1);
                int left_seq2 = (c2_left == 1) ? static_cast<int>(x2) - 1 : static_cast<int>(x2);

                // D is indexed by the offset of that position from the arc endpoint
                int x1_offset = left_seq1 - static_cast<int>(al);
                int x2_offset = left_seq2 - static_cast<int>(bl);

                // Check offset bounds
                if (std::abs(x1_offset) > static_cast<int>(max_shifts_) ||
                    std::abs(x2_offset) > static_cast<int>(max_shifts_))
                    continue;

                // Compute left column score (sequence score + shift penalty)
                score_t left_seq_score = score_t(0);
                if (c1_left == 1 && c2_left == 1)
                    left_seq_score = score_t(locarna_scoring_->basematch(left_seq1, left_seq2));
                else if (c1_left == 1)
                    left_seq_score = score_t(locarna_scoring_->gapA(left_seq1));
                else if (c2_left == 1)
                    left_seq_score = score_t(locarna_scoring_->gapB(left_seq2));
                // else both gaps: score remains 0

                // Shift penalty: structure always has (1,1) at arc endpoints
                score_t left_shift = shift_scoring_.shift_penalty(c1_left, c2_left, 1, 1);
                score_t left_total = left_seq_score + left_shift;

                // Iterate over right endpoint shifts
                for (int shift_y1 = -static_cast<int>(max_shifts_);
                     shift_y1 <= static_cast<int>(max_shifts_); ++shift_y1) {
                    for (int shift_y2 = -static_cast<int>(max_shifts_);
                         shift_y2 <= static_cast<int>(max_shifts_); ++shift_y2) {

                        int y1 = static_cast<int>(ar) + shift_y1;
                        int y2 = static_cast<int>(br) + shift_y2;

                        // Bounds checks for right sequence positions
                        if (y1 <= static_cast<int>(x1) || y2 <= static_cast<int>(x2))
                            continue;
                        if (y1 > static_cast<int>(lenA) || y2 > static_cast<int>(lenB))
                            continue;

                        // These shifts become the D matrix offsets at right end
                        int y1_offset = shift_y1;
                        int y2_offset = shift_y2;

                        // Iterate over right column gap patterns (4 patterns)
                        for (int c1_right = 0; c1_right <= 1; ++c1_right) {
                            for (int c2_right = 0; c2_right <= 1; ++c2_right) {

                                // Determine M matrix predecessor based on right gap pattern
                                // If c1_right=1, sequence position y1 is consumed, so M predecessor at y1-1
                                // If c1_right=0, gap in first sequence, M predecessor at y1
                                int m_y1 = (c1_right == 1) ? y1 - 1 : y1;
                                int m_y2 = (c2_right == 1) ? y2 - 1 : y2;

                                // Bounds check for M access
                                // Lower bound: must be >= x1-1, x2-1 (base case position)
                                if (m_y1 < static_cast<int>(x1) - 1 || m_y2 < static_cast<int>(x2) - 1)
                                    continue;

                                // Upper bound: M matrix interior ends at (ar-1, br-1), so max sequence
                                // positions filled are (ar-1) + max_shifts and (br-1) + max_shifts
                                int max_m_y1 = static_cast<int>(ar_interior) + static_cast<int>(max_shifts_);
                                int max_m_y2 = static_cast<int>(br_interior) + static_cast<int>(max_shifts_);
                                if (m_y1 > max_m_y1 || m_y2 > max_m_y2)
                                    continue;

                                // Get M score at the right interior boundary
                                score_t m_score = M_->get(m_y1, m_y2, ar_interior, br_interior);
                                if (m_score == score_t::neg_infty)
                                    continue;

                                // Compute right column score (sequence score + shift penalty)
                                score_t right_seq_score = score_t(0);
                                if (c1_right == 1 && c2_right == 1)
                                    right_seq_score = score_t(locarna_scoring_->basematch(y1, y2));
                                else if (c1_right == 1)
                                    right_seq_score = score_t(locarna_scoring_->gapA(y1));
                                else if (c2_right == 1)
                                    right_seq_score = score_t(locarna_scoring_->gapB(y2));
                                // else both gaps: score remains 0

                                score_t right_shift = shift_scoring_.shift_penalty(c1_right, c2_right, 1, 1);
                                score_t right_total = right_seq_score + right_shift;

                                // Compute total score for this configuration
                                score_t total = left_total + m_score + right_total + arc_score;

                                // Convert offsets back to sequence positions for D matrix access
                                // D matrix set() method expects sequence positions, not offsets
                                int x1_seq = static_cast<int>(al) + x1_offset;
                                int x2_seq = static_cast<int>(bl) + x2_offset;
                                int y1_seq = static_cast<int>(ar) + y1_offset;
                                int y2_seq = static_cast<int>(br) + y2_offset;

                                // Get existing D entry (if any)
                                score_t existing = D_->get(arcA, arcB, x1_seq, x2_seq, y1_seq, y2_seq);

                                // Store if this is better (or if entry doesn't exist)
                                if (existing == score_t::neg_infty || total > existing) {
                                    D_->set(arcA, arcB, x1_seq, x2_seq, y1_seq, y2_seq, total);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



ShiftAligner::score_t
ShiftAligner::compute_unpaired_score(size_type x1, size_type x2,
                                      size_type x3, size_type x4,
                                      size_type y1, size_type y2,
                                      size_type y3, size_type y4) {
    // Case 1: the column at (y1,y2,y3,y4) is unpaired. Try each of the 15 valid
    // column types c, score its predecessor y-c, and keep the best:
    //
    //     M(y) = max over c of  M(y - c) + u_s + v_s + w_s
    //
    // u_s scores the sequence layer, v_s the structure layer (gaps only — a
    // structure match scores 0 here; base pairs are rewarded via the arc-match
    // score in Case 2), and w_s penalizes a disagreement between the two layers.

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
        // Base case is at (x1-1, x2-1, x3-1, x4-1), so allow access to those positions
        if (px1 < static_cast<int>(x1)-1 || px2 < static_cast<int>(x2)-1 ||
            px3 < static_cast<int>(x3)-1 || px4 < static_cast<int>(x4)-1)
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
        if (col.c1 == 1 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->basematch(y1, y2));
        } else if (col.c1 == 1 && col.c2 == 0) {
            u_s = score_t(locarna_scoring_->gapA(y1));
        } else if (col.c1 == 0 && col.c2 == 1) {
            u_s = score_t(locarna_scoring_->gapB(y2));
        }
        // else: (0,0) gap-on-gap, u_s = 0

        // v_s: structure score for V layer (positions y3, y4)
        score_t v_s(0);
        if (col.c3 == 1 && col.c4 == 1) {
            // v_s = score_t(locarna_scoring_->basematch(y3, y4)); // 0 for match/mismatch
        } else if (col.c3 == 1 && col.c4 == 0) {
            v_s = score_t(locarna_scoring_->gapA(y3));
        } else if (col.c3 == 0 && col.c4 == 1) {
            v_s = score_t(locarna_scoring_->gapB(y4));
        }
        // else: (0,0) gap-on-gap, v_s = 0

        // w_s: shift penalty using column components directly
        score_t w_s = shift_scoring_.shift_penalty(col.c1, col.c2, col.c3, col.c4);

        score_t total_score = pred_score + u_s + v_s + w_s;

        if (total_score > best_score) {
            best_score = total_score;
        }
    }

    return best_score;
}

ShiftAligner::score_t
ShiftAligner::compute_paired_score(size_type x1, size_type x2,
                                    size_type x3, size_type x4,
                                    size_type y1, size_type y2,
                                    size_type y3, size_type y4) {
    // Case 2: the structure layer closes an arc match at (y3, y4). For each such
    // arc, try every sequence position (z1, z2) its left ends could sit at, and
    // combine the alignment before the arc with the arc's complete score:
    //
    //     M(y) = max over (arc, z) of  M(z1-1, z2-1, z3-1, z4-1) + D(arc, z, y)
    //
    // D already bundles both endpoint columns, the interior, and the arc-match
    // reward.

    score_t best_score = score_t::neg_infty;

    // Get all arc matches with right ends at current structure positions (y3, y4)
    const LocARNA::ArcMatchIdxVec& arcmatches_at_right =
        arc_matches_->common_right_end_list(y3, y4);

    for (LocARNA::ArcMatch::idx_type am_idx : arcmatches_at_right) {
        const LocARNA::ArcMatch& arc_match = arc_matches_->arcmatch(am_idx);
        const LocARNA::BasePairs__Arc& arcA = arc_match.arcA();
        const LocARNA::BasePairs__Arc& arcB = arc_match.arcB();

        size_type z3 = arcA.left();   // Left arc endpoint in structure A
        size_type z4 = arcB.left();   // Left arc endpoint in structure B

        // Arc must be inside current region (structure boundaries are at x3-1, x4-1)
        if (z3 < x3 || z4 < x4)
            continue;

        // Iterate over all possible sequence positions z1, z2 at arc left ends
        for (int shift_z1 = -static_cast<int>(max_shifts_);
             shift_z1 <= static_cast<int>(max_shifts_); ++shift_z1) {
            for (int shift_z2 = -static_cast<int>(max_shifts_);
                 shift_z2 <= static_cast<int>(max_shifts_); ++shift_z2) {

                int z1 = static_cast<int>(z3) + shift_z1;
                int z2 = static_cast<int>(z4) + shift_z2;

                // Check bounds for z positions (must be within region)
                if (z1 < static_cast<int>(x1) || z2 < static_cast<int>(x2))
                    continue;

                // Predecessor position: just before arc left endpoints
                int px1 = z1 - 1;  // Sequence advances from px1 to z1
                int px2 = z2 - 1;  // Sequence advances from px2 to z2
                int px3 = static_cast<int>(z3) - 1;  // Structure advances from px3 to z3
                int px4 = static_cast<int>(z4) - 1;  // Structure advances from px4 to z4

                // Check bounds (base case is at x1-1, x2-1, x3-1, x4-1)
                if (px1 < static_cast<int>(x1)-1 || px2 < static_cast<int>(x2)-1 ||
                    px3 < static_cast<int>(x3)-1 || px4 < static_cast<int>(x4)-1)
                    continue;

                // Check δ_max constraint at predecessor
                if (std::abs(px1 - px3) > static_cast<int>(max_shifts_))
                    continue;
                if (std::abs(px2 - px4) > static_cast<int>(max_shifts_))
                    continue;

                // Get M score before the arc
                score_t m_before = M_->get(px1, px2, px3, px4);

                // Check if D matrix access would be valid
                // D matrix offsets must be within [-max_shifts, +max_shifts]
                // Offsets are computed relative to arc endpoints
                int d_x_off1 = z1 - static_cast<int>(arcA.left());
                int d_x_off2 = z2 - static_cast<int>(arcB.left());
                int d_y_off1 = static_cast<int>(y1) - static_cast<int>(arcA.right());
                int d_y_off2 = static_cast<int>(y2) - static_cast<int>(arcB.right());

                // Skip if the D entry would lie outside the stored offset range
                if (std::abs(d_x_off1) > static_cast<int>(max_shifts_) ||
                    std::abs(d_x_off2) > static_cast<int>(max_shifts_) ||
                    std::abs(d_y_off1) > static_cast<int>(max_shifts_) ||
                    std::abs(d_y_off2) > static_cast<int>(max_shifts_))
                    continue;

                // Check if offset matrix is allocated for this arc pair
                if (!D_->is_offsetmatrix_allocated(arcA.idx(), arcB.idx()))
                    continue;

                // D score includes everything:
                // - Interior M score
                // - Arc match score (structure probabilities)
                // - Left endpoint column (sequence score + shift penalty)
                // - Right endpoint column (sequence score + shift penalty)
                score_t d_score = D_->get(arcA, arcB, z1, z2, y1, y2);

                // Total score: just combine predecessor with complete arc score
                score_t total_score = m_before + d_score;

                if (total_score > best_score) {
                    best_score = total_score;
                }
            }
        }
    }

    return best_score;
}

} // namespace RNAShiftAlign
