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

TEST_CASE("ShiftAligner - Stub methods", "[shift_aligner]") {
    ShiftAligner aligner("GCGC", "GCGC");

    SECTION("align() returns a score (currently stub)") {
        // Currently fill_M_unpaired is empty, so this will just
        // check that the method runs without crashing
        score_t score = aligner.align();
        // Empty implementation should return the base case M(0,0,0,0) = 0
        // But M(len,len,len,len) is uninitialized (-inf), so expect that
        CHECK(score.is_neg_infty());
    }

    SECTION("traceback() runs without crashing (currently stub)") {
        // Should not throw even though it's a stub
        REQUIRE_NOTHROW(aligner.traceback());
    }

    SECTION("get_score() before align() returns initial value") {
        CHECK(aligner.get_score() == score_t(0));
    }
}
