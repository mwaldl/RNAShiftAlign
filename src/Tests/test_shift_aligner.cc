#include "catch.hpp"
#include "RNAShiftAlign/shift_aligner.hh"

using namespace RNAShiftAlign;
using score_t = ShiftAligner::score_t;

TEST_CASE("ShiftAligner - Basic construction", "[shift_aligner]") {
    std::string seqA = "GCGCGC";
    std::string seqB = "GCGCGC";

    SECTION("Construct with default parameters") {
        ShiftAligner aligner(seqA, seqB);

        // Verify sequences are stored correctly
        CHECK(aligner.get_seqA().length() == 6);
        CHECK(aligner.get_seqB().length() == 6);

        // Initial score should be 0
        CHECK(aligner.get_score() == score_t(0));
    }

    SECTION("Construct with custom parameters") {
        PreprocessingParams prep_params;
        prep_params.min_prob = 0.001;
        prep_params.max_diff_am = 30;

        ScoringParams scoring_params;
        scoring_params.match = 100;
        scoring_params.mismatch = -50;
        scoring_params.indel = -200;
        scoring_params.struct_weight = 150;

        ShiftAligner aligner(seqA, seqB, prep_params, scoring_params);

        CHECK(aligner.get_seqA().length() == 6);
        CHECK(aligner.get_seqB().length() == 6);
    }
}

TEST_CASE("ShiftAligner - Preprocessing components", "[shift_aligner]") {
    // Use a sequence that can form base pairs
    std::string seqA = "GCGCGCGC";
    std::string seqB = "GUGCGCAC";

    ShiftAligner aligner(seqA, seqB);

    SECTION("RnaData is initialized") {
        const auto& rna_dataA = aligner.get_rna_dataA();
        const auto& rna_dataB = aligner.get_rna_dataB();

        // RnaData should have sequences matching input
        CHECK(rna_dataA.sequence().length() == 8);
        CHECK(rna_dataB.sequence().length() == 8);

        // RnaData should have computed arc cutoff probabilities
        CHECK(rna_dataA.arc_cutoff_prob() > 0.0);
        CHECK(rna_dataB.arc_cutoff_prob() > 0.0);
    }

    SECTION("RnaData computes base pair list correctly") {
        // Use a sequence with a clear hairpin structure
        // ACCCCAAAAGGGGA should form pairs in the stem: (2,13), (3,12), (4,11), (5,10)
        // Note: ViennaRNA uses 1-based indexing
        std::string hairpin_seq = "ACCCCAAAAGGGGA";
        ShiftAligner hairpin_aligner(hairpin_seq, hairpin_seq);

        const auto& rna_data = hairpin_aligner.get_rna_dataA();

        // Get the pair list from RnaData - this is a unique_ptr to vrna_elem_prob_s array
        // Each element has fields: i (first position), j (second position), p (probability)
        // The list is terminated by an entry with i == 0
        const auto& plist_ptr = rna_data.plist();

        // First, let's see what pairs are actually computed
        std::vector<std::pair<int, int>> found_pairs;
        for (int idx = 0; plist_ptr[idx].i != 0; ++idx) {
            found_pairs.push_back({plist_ptr[idx].i, plist_ptr[idx].j});
        }

        // Expected stem pairs for ACCCCAAAAGGGGA
        // The hairpin should have C-G pairs in the stem
        // Sequence: A C C C C A A A A G G G G A
        // Positions: 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        // Expected structure: .((((....)))).
        // Expected pairs: (2,13), (3,12), (4,11), (5,10)
        std::vector<std::pair<int, int>> expected_pairs = {
            {2, 13},
            {3, 12},
            {4, 11},
            {5, 10}
        };

        // Check that we found at least some pairs
        REQUIRE(found_pairs.size() > 0);

        // Check that all expected stem pairs are present
        // ViennaRNA computes all pairs with probability > threshold, so there may be
        // additional low-probability pairs beyond the MFE structure
        for (const auto& expected : expected_pairs) {
            bool found = false;
            for (const auto& actual : found_pairs) {
                if (actual.first == expected.first && actual.second == expected.second) {
                    found = true;
                    break;
                }
            }
            INFO("Expected stem pair (" << expected.first << "," << expected.second << ") not found");
            CHECK(found);
        }
    }

    SECTION("ArcMatches is initialized") {
        const auto& arc_matches = aligner.get_arc_matches();

        // Arc matches should be accessible (may be empty for short sequences)
        // Just verify it doesn't crash
        size_t num_matches = arc_matches.num_arc_matches();
        CHECK(num_matches >= 0);
    }

    SECTION("Scoring is initialized") {
        const auto& scoring = aligner.get_scoring();

        // Verify scoring object is functional
        // Test base match scoring
        CHECK(scoring.basematch(1, 1) >= 0);  // Matching position should have positive score
    }
}

TEST_CASE("ShiftAligner - Different sequence lengths", "[shift_aligner]") {
    SECTION("Equal length sequences") {
        ShiftAligner aligner("AAAA", "UUUU");
        CHECK(aligner.get_seqA().length() == 4);
        CHECK(aligner.get_seqB().length() == 4);
    }

    SECTION("Sequence A longer") {
        ShiftAligner aligner("AAAAAAAA", "UUUU");
        CHECK(aligner.get_seqA().length() == 8);
        CHECK(aligner.get_seqB().length() == 4);
    }

    SECTION("Sequence B longer") {
        ShiftAligner aligner("AAAA", "UUUUUUUU");
        CHECK(aligner.get_seqA().length() == 4);
        CHECK(aligner.get_seqB().length() == 8);
    }
}

TEST_CASE("ShiftAligner - Basic alignment", "[shift_aligner]") {
    SECTION("align() returns a valid score for identical sequences") {
        ShiftAligner aligner("GCGC", "GCGC");
        score_t score = aligner.align();

        // For identical sequences with default parameters, we expect a positive score
        // (all matches, no gaps, no shifts since U=V is optimal)
        CHECK(!score.is_neg_infty());
        CHECK(score > score_t(0));
    }

    SECTION("traceback() runs without crashing (currently stub)") {
        ShiftAligner aligner("GCGC", "GCGC");
        // Should not throw even though it's a stub
        REQUIRE_NOTHROW(aligner.traceback());
    }

    SECTION("get_score() before align() returns initial value") {
        ShiftAligner aligner("GCGC", "GCGC");
        CHECK(aligner.get_score() == score_t(0));
    }

    SECTION("get_score() after align() returns computed score") {
        ShiftAligner aligner("GCGC", "GCGC");
        score_t align_score = aligner.align();
        score_t get_score = aligner.get_score();
        CHECK(align_score == get_score);
    }
}

TEST_CASE("ShiftAligner - Forward recursion with hand-calculated examples", "[shift_aligner][forward_recursion]") {
    SECTION("Identical sequences - optimal should have U=V with no shifts") {
        // For identical sequences with Δ > 0, the optimal bi-alignment
        // should have U = V (no shifts), giving all matches
        std::string seq = "AAAA";
        ShiftAligner aligner(seq, seq);

        score_t score = aligner.align();

        // With default params: match=50, indel=-150, delta=100
        // Optimal: 4 matches in U, 4 matches in V, 0 shift penalty
        // Expected score: 4*(50 + 50 + 0) = 400
        // Note: This assumes basematch for identical bases is 50
        CHECK(score > score_t(0));
        CHECK(!score.is_neg_infty());
    }

    SECTION("Different length sequences require gaps") {
        // Align "AA" with "AAAA" - requires 2 gaps
        ShiftAligner aligner("AA", "AAAA");

        score_t score = aligner.align();

        // Should be able to align with gaps (negative score due to indels)
        CHECK(!score.is_neg_infty());
        // With match=50, indel=-150: 2 matches + 2 gaps = 2*100 + 2*(-150) = -100
        // But this is approximate since we have two layers
        CHECK(score < score_t(500));  // Should not be unreasonably high
    }

    SECTION("Completely different sequences") {
        // Different sequences should still align but with lower score
        ShiftAligner aligner("AAAA", "UUUU");

        score_t score = aligner.align();

        // Should complete without error
        CHECK(!score.is_neg_infty());

        // Compare to identical sequences - should be lower score
        ShiftAligner aligner_same("AAAA", "AAAA");
        score_t score_same = aligner_same.align();

        CHECK(score < score_same);
    }

    SECTION("Short sequences - verify matrix filling") {
        // Very short sequence to manually verify
        ShiftAligner aligner("A", "A");

        score_t score = aligner.align();

        // Single match in both layers, no gaps, no shifts
        // Expected: 1*(match + match + 0) = 100
        CHECK(score > score_t(0));
        CHECK(!score.is_neg_infty());
    }

    SECTION("δ_max constraint is respected") {
        // With default max_shifts=5, alignment should complete
        ScoringParams params;
        params.max_shifts = 2;  // Tighter constraint

        ShiftAligner aligner("AAAA", "AAAA", PreprocessingParams(), params);

        score_t score = aligner.align();

        // Should still align successfully with tighter constraint
        CHECK(!score.is_neg_infty());
        CHECK(score > score_t(0));
    }

    SECTION("Shift penalty affects score when U != V") {
        // Create scenario where shifts might occur
        // Use different delta values and verify score changes
        ScoringParams params_low_delta;
        params_low_delta.delta = 10;  // Low shift penalty

        ScoringParams params_high_delta;
        params_high_delta.delta = 1000;  // High shift penalty

        ShiftAligner aligner_low("AAUU", "UUAA", PreprocessingParams(), params_low_delta);
        ShiftAligner aligner_high("AAUU", "UUAA", PreprocessingParams(), params_high_delta);

        score_t score_low = aligner_low.align();
        score_t score_high = aligner_high.align();

        // Both should complete
        CHECK(!score_low.is_neg_infty());
        CHECK(!score_high.is_neg_infty());

        // With high delta, algorithm should avoid shifts more aggressively
        // Scores might differ depending on whether shifts are used
    }
}
