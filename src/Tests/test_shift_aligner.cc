#include "catch.hpp"
#include "RNAShiftAlign/shift_aligner.hh"

using namespace RNAShiftAlign;
using score_t = ShiftAligner::score_t;

TEST_CASE("ShiftAligner - Constructor and initialization", "[shift_aligner]") {
    std::string seqA = "GCGCGC";
    std::string seqB = "GCGCGC";

    SECTION("Construct with simple sequences") {
        ShiftAligner aligner(seqA, seqB, score_t(100), 2);

        CHECK(aligner.get_seqA().length() == 6);
        CHECK(aligner.get_seqB().length() == 6);
        CHECK(aligner.get_delta() == score_t(100));
        CHECK(aligner.get_max_shifts() == 2);
    }

    SECTION("Construct with different max_shifts") {
        ShiftAligner aligner1(seqA, seqB, score_t(100), 0);
        ShiftAligner aligner2(seqA, seqB, score_t(100), 5);

        CHECK(aligner1.get_max_shifts() == 0);
        CHECK(aligner2.get_max_shifts() == 5);
    }

    SECTION("Construct with different delta values") {
        ShiftAligner aligner1(seqA, seqB, score_t(50), 2);
        ShiftAligner aligner2(seqA, seqB, score_t(200), 2);

        CHECK(aligner1.get_delta() == score_t(50));
        CHECK(aligner2.get_delta() == score_t(200));
    }
}

TEST_CASE("ShiftAligner - Different sequence lengths", "[shift_aligner]") {
    SECTION("Equal length sequences") {
        ShiftAligner aligner("AAAA", "UUUU", score_t(100), 2);
        CHECK(aligner.get_seqA().length() == 4);
        CHECK(aligner.get_seqB().length() == 4);
    }

    SECTION("Sequence A longer") {
        ShiftAligner aligner("AAAAAAAA", "UUUU", score_t(100), 2);
        CHECK(aligner.get_seqA().length() == 8);
        CHECK(aligner.get_seqB().length() == 4);
    }

    SECTION("Sequence B longer") {
        ShiftAligner aligner("AAAA", "UUUUUUUU", score_t(100), 2);
        CHECK(aligner.get_seqA().length() == 4);
        CHECK(aligner.get_seqB().length() == 8);
    }
}

TEST_CASE("ShiftAligner - Sequence content", "[shift_aligner]") {

    SECTION("Sequences are correctly stored") {
        std::string seqA = "GCAUCG";
        std::string seqB = "GUAUCG";

        ShiftAligner aligner(seqA, seqB, score_t(100), 2);

        // Sequences should be stored with correct lengths
        CHECK(aligner.get_seqA().length() == 6);
        CHECK(aligner.get_seqB().length() == 6);

        // Verify we can access the sequence objects
        const auto& seq_a = aligner.get_seqA();
        const auto& seq_b = aligner.get_seqB();

        // Both sequences should be valid
        CHECK(seq_a.length() > 0);
        CHECK(seq_b.length() > 0);
    }
}

TEST_CASE("ShiftAligner - Stub methods", "[shift_aligner]") {
    ShiftAligner aligner("GCGC", "GCGC", score_t(100), 2);

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

TEST_CASE("ShiftAligner - Zero max_shifts edge case", "[shift_aligner]") {

    SECTION("max_shifts = 0 means no shifts allowed") {
        ShiftAligner aligner("GC", "GC", score_t(100), 0);
        CHECK(aligner.get_max_shifts() == 0);

        // With max_shifts = 0, only M(i,j,i,j) entries are valid
        // This means U and V must be identical (no shift allowed)
    }
}
