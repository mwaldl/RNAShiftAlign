#include "catch.hpp"

#include <ShiftAlign/shift_aligner.hh>
#include <ShiftAlign/shift_scoring.hh>

using namespace RNAShiftAlign;
using score_t = ShiftAligner::score_t;

/**
 * @file Tests for structure-aware shift alignment
 *
 * Tests the full algorithm including:
 * - D-matrix computation (fill_D)
 * - Local M-matrix filling (fill_M_local)
 * - Case 2 recursion (paired positions)
 * - Integration of structure scoring
 */

TEST_CASE("Simple hairpin - identical sequences", "[structure][hairpin]") {
    // Test Case 8: Simple hairpin end-to-end
    //
    // Sequence: CCCAAAGGG
    // Structure:
    //   CCC AAA GGG
    //   (((  )  )))
    //   123     987
    //
    // Expected base pairs (1-indexed):
    // - (1,9): C-G
    // - (2,8): C-G
    // - (3,7): C-G
    //
    // This forms a stem with 3 base pairs and a 3-nt loop (AAA)

    std::string seq = "CCCAAAGGG";

    SECTION("Structure alignment produces higher score than unpaired") {
        PreprocessingParams prep_struct;
        prep_struct.min_prob = 0.0001;  // Low threshold to catch all base pairs

        PreprocessingParams prep_no_struct;
        prep_no_struct.min_prob = 0.99;  // High threshold - no base pairs accepted

        ScoringParams scoring;
        scoring.match = 50;
        scoring.mismatch = 0;
        scoring.indel = -100;
        scoring.indel_opening = -500;
        scoring.struct_weight = 200;  // Strong structure bonus
        scoring.tau_factor = 0;
        scoring.delta = -150;  // Shift penalty
        scoring.max_shifts = 2;

        // Align with structure
        ShiftAligner aligner_struct(seq, seq, prep_struct, scoring);
        auto score_struct = aligner_struct.align();

        // Align without structure (unpaired only)
        ShiftAligner aligner_no_struct(seq, seq, prep_no_struct, scoring);
        auto score_no_struct = aligner_no_struct.align();

        // Structure version should have higher score due to arc match bonuses
        INFO("Score with structure: " << score_struct);
        INFO("Score without structure: " << score_no_struct);
        INFO("Number of arc matches: " << aligner_struct.get_arc_matches().num_arc_matches());

        REQUIRE(aligner_struct.get_arc_matches().num_arc_matches() > 0);
        REQUIRE(score_struct > score_no_struct);
    }

    SECTION("Verify arc matches are found") {
        PreprocessingParams prep;
        prep.min_prob = 0.0001;

        ScoringParams scoring;
        scoring.match = 50;
        scoring.indel = -100;
        scoring.struct_weight = 200;
        scoring.delta = -150;
        scoring.max_shifts = 2;

        ShiftAligner aligner(seq, seq, prep, scoring);

        const auto& arc_matches = aligner.get_arc_matches();
        size_t num_arcs = arc_matches.num_arc_matches();

        INFO("Found " << num_arcs << " arc matches");

        // Should find at least the 3 stem base pairs
        // (Might find more due to stacking or alternative structures)
        REQUIRE(num_arcs >= 3);

        // Check that some arcs are in the expected positions
        bool found_outer = false;  // (1,9)
        bool found_middle = false; // (2,8)
        bool found_inner = false;  // (3,7)

        for (size_t i = 0; i < num_arcs; ++i) {
            const auto& am = arc_matches.arcmatch(i);
            auto arcA = am.arcA();
            auto arcB = am.arcB();

            INFO("Arc match " << i << ": (" << arcA.left() << "," << arcA.right()
                 << ") x (" << arcB.left() << "," << arcB.right() << ")");

            // Check if this matches our expected arcs
            if (arcA.left() == 1 && arcA.right() == 9 &&
                arcB.left() == 1 && arcB.right() == 9) {
                found_outer = true;
            }
            if (arcA.left() == 2 && arcA.right() == 8 &&
                arcB.left() == 2 && arcB.right() == 8) {
                found_middle = true;
            }
            if (arcA.left() == 3 && arcA.right() == 7 &&
                arcB.left() == 3 && arcB.right() == 7) {
                found_inner = true;
            }
        }

        // Should find at least one of the expected arcs
        REQUIRE((found_outer || found_middle || found_inner));
    }

    SECTION("Alignment score is reasonable") {
        PreprocessingParams prep;
        prep.min_prob = 0.0001;

        ScoringParams scoring;
        scoring.match = 50;
        scoring.mismatch = 0;
        scoring.indel = -100;
        scoring.struct_weight = 200;
        scoring.delta = -150;
        scoring.max_shifts = 2;

        ShiftAligner aligner(seq, seq, prep, scoring);
        auto score = aligner.align();

        // Expected score breakdown:
        // - 9 base matches: 9 * 50 = 450
        // - 3 arc matches: 3 * struct_weight ≈ 600 (depends on probabilities)
        // - No shifts needed (identical sequences): 0
        // - No gaps: 0
        //
        // Minimum expected: 450 (just base matches)
        // Maximum expected: ~1050 (base matches + arc bonuses)

        INFO("Final score: " << score);

        REQUIRE(score >= 450);  // At least base matches
        REQUIRE(score <= 2000); // Sanity check - not wildly wrong
    }
}

TEST_CASE("D-matrix extraction", "[structure][fill_D]") {
    // Test Case 5: Verify D-matrix is filled correctly
    //
    // Use same sequence: CCCAAAGGG
    // We can't directly inspect D-matrix values without adding accessors,
    // but we can verify the algorithm completes and produces reasonable results

    std::string seq = "CCCAAAGGG";

    PreprocessingParams prep;
    prep.min_prob = 0.0001;

    ScoringParams scoring;
    scoring.match = 50;
    scoring.indel = -100;
    scoring.struct_weight = 200;
    scoring.delta = -150;
    scoring.max_shifts = 2;

    ShiftAligner aligner(seq, seq, prep, scoring);

    SECTION("Algorithm completes without errors") {
        // This calls fill_D() internally
        REQUIRE_NOTHROW(aligner.align());
    }

    SECTION("Score indicates D-matrix was used") {
        auto score = aligner.align();

        // Score should be higher than just base matches (450)
        // indicating arc matches from D-matrix were used
        REQUIRE(score > 450);
    }
}

TEST_CASE("Nested arcs", "[structure][nested]") {
    // Test Case 6: Verify inner arcs are processed before outer
    //
    // Sequence: CCCAAAGGG has nested structure
    // - Outer: (1,9)
    // - Middle: (2,8)
    // - Inner: (3,7)
    //
    // When filling D for outer arc, D for inner arcs should already exist

    std::string seq = "CCCAAAGGG";

    PreprocessingParams prep;
    prep.min_prob = 0.0001;

    ScoringParams scoring;
    scoring.match = 50;
    scoring.indel = -100;
    scoring.struct_weight = 200;
    scoring.delta = -150;
    scoring.max_shifts = 2;

    ShiftAligner aligner(seq, seq, prep, scoring);

    SECTION("Nested structure is handled") {
        // Algorithm should handle nested arcs correctly
        REQUIRE_NOTHROW(aligner.align());

        auto score = aligner.align();

        // With nested arcs, score should reflect:
        // - Base matches in loop (AAA): 3 * 50 = 150
        // - Base matches in stem (CCC, GGG): 6 * 50 = 300
        // - Arc match bonuses for each base pair
        //
        // Total should be substantially higher than unpaired
        REQUIRE(score > 450);
    }
}

TEST_CASE("No arc matches - fallback to unpaired", "[structure][edge_case]") {
    // Test Case 11: Sequences with no common arc matches
    //
    // Different sequences that won't form matching base pairs

    std::string seqA = "AAAAAAA";  // All A - no structure
    std::string seqB = "CCCCCCC";  // All C - no structure

    PreprocessingParams prep;
    prep.min_prob = 0.0001;

    ScoringParams scoring;
    scoring.match = 50;
    scoring.mismatch = 0;
    scoring.indel = -100;
    scoring.struct_weight = 200;
    scoring.delta = -150;
    scoring.max_shifts = 2;

    ShiftAligner aligner(seqA, seqB, prep, scoring);

    SECTION("No arc matches found") {
        const auto& arc_matches = aligner.get_arc_matches();
        INFO("Number of arc matches: " << arc_matches.num_arc_matches());

        // Should have 0 or very few arc matches
        REQUIRE(arc_matches.num_arc_matches() == 0);
    }

    SECTION("Algorithm handles empty D-matrix gracefully") {
        // Should complete without errors even with no arc matches
        REQUIRE_NOTHROW(aligner.align());

        // Score should only reflect base mismatches (all 0) and gaps
        // or mismatches depending on how sequences are aligned
        auto score = aligner.align();
        INFO("Score with no arc matches: " << score);

        // Score should be 0 or negative (no matches between A and C)
        REQUIRE(score <= 0);
    }
}

TEST_CASE("Delta max constraint with structure", "[structure][delta_max]") {
    // Test Case 12: Verify shift constraints work with arc matches

    std::string seq = "CCCAAAGGG";

    PreprocessingParams prep;
    prep.min_prob = 0.0001;

    ScoringParams scoring;
    scoring.match = 50;
    scoring.indel = -100;
    scoring.struct_weight = 200;
    scoring.delta = -150;

    SECTION("Very tight delta_max constraint") {
        scoring.max_shifts = 0;  // No shifts allowed

        ShiftAligner aligner(seq, seq, prep, scoring);
        auto score = aligner.align();

        // With identical sequences and no shifts allowed,
        // should get perfect alignment (all matches + arcs)
        REQUIRE(score > 450);
    }

    SECTION("Looser delta_max constraint") {
        scoring.max_shifts = 3;  // More shifts allowed

        ShiftAligner aligner(seq, seq, prep, scoring);
        auto score = aligner.align();

        // Score should be similar (identical sequences don't need shifts)
        REQUIRE(score > 450);
    }
}
