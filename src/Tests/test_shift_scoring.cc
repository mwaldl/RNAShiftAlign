#include "catch.hpp"
#include "ShiftAlign/shift_scoring.hh"

using namespace RNAShiftAlign;
using score_t = ShiftScoring::score_t;

TEST_CASE("ShiftScoring - Constructor and basic properties", "[shift_scoring]") {
    SECTION("Constructor accepts delta parameter") {
        ShiftScoring scorer(score_t(100));
        CHECK(scorer.get_delta() == score_t(100));
    }

    SECTION("Different delta values") {
        ShiftScoring scorer1(score_t(50));
        ShiftScoring scorer2(score_t(200));
        CHECK(scorer1.get_delta() == score_t(50));
        CHECK(scorer2.get_delta() == score_t(200));
    }
}

TEST_CASE("ShiftScoring - Diagonal (identical column types)", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    SECTION("MATCH vs MATCH") {
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::MATCH) == score_t(0));
    }

    SECTION("DEL_A vs DEL_A") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::DEL_A) == score_t(0));
    }

    SECTION("INS_A vs INS_A") {
        CHECK(scorer.shift_penalty(ColumnType::INS_A, ColumnType::INS_A) == score_t(0));
    }

    SECTION("DEL_B vs DEL_B") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_B, ColumnType::DEL_B) == score_t(0));
    }

    SECTION("INS_B vs INS_B") {
        CHECK(scorer.shift_penalty(ColumnType::INS_B, ColumnType::INS_B) == score_t(0));
    }
}

TEST_CASE("ShiftScoring - MATCH vs gaps (penalty Δ)", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    SECTION("MATCH in U, DEL_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A) == score_t(100));
    }

    SECTION("MATCH in U, INS_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::INS_A) == score_t(100));
    }

    SECTION("MATCH in U, DEL_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_B) == score_t(100));
    }

    SECTION("MATCH in U, INS_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::INS_B) == score_t(100));
    }

    SECTION("DEL_A in U, MATCH in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::MATCH) == score_t(100));
    }

    SECTION("INS_A in U, MATCH in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_A, ColumnType::MATCH) == score_t(100));
    }

    SECTION("DEL_B in U, MATCH in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_B, ColumnType::MATCH) == score_t(100));
    }

    SECTION("INS_B in U, MATCH in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_B, ColumnType::MATCH) == score_t(100));
    }
}

TEST_CASE("ShiftScoring - Same gap types (penalty 0)", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    SECTION("DEL_A in U, DEL_B in V (both deletions)") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::DEL_B) == score_t(0));
    }

    SECTION("DEL_B in U, DEL_A in V (both deletions)") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_B, ColumnType::DEL_A) == score_t(0));
    }

    SECTION("INS_A in U, INS_B in V (both insertions)") {
        CHECK(scorer.shift_penalty(ColumnType::INS_A, ColumnType::INS_B) == score_t(0));
    }

    SECTION("INS_B in U, INS_A in V (both insertions)") {
        CHECK(scorer.shift_penalty(ColumnType::INS_B, ColumnType::INS_A) == score_t(0));
    }
}

TEST_CASE("ShiftScoring - Different gap types (penalty 2Δ)", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    SECTION("DEL_A in U, INS_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_A) == score_t(200));
    }

    SECTION("DEL_A in U, INS_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_B) == score_t(200));
    }

    SECTION("DEL_B in U, INS_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_B, ColumnType::INS_A) == score_t(200));
    }

    SECTION("DEL_B in U, INS_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::DEL_B, ColumnType::INS_B) == score_t(200));
    }

    SECTION("INS_A in U, DEL_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_A, ColumnType::DEL_A) == score_t(200));
    }

    SECTION("INS_A in U, DEL_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_A, ColumnType::DEL_B) == score_t(200));
    }

    SECTION("INS_B in U, DEL_A in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_B, ColumnType::DEL_A) == score_t(200));
    }

    SECTION("INS_B in U, DEL_B in V") {
        CHECK(scorer.shift_penalty(ColumnType::INS_B, ColumnType::DEL_B) == score_t(200));
    }
}

TEST_CASE("ShiftScoring - Symmetry", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    SECTION("Penalty is symmetric for MATCH vs gap") {
        auto p1 = scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A);
        auto p2 = scorer.shift_penalty(ColumnType::DEL_A, ColumnType::MATCH);
        CHECK(p1 == p2);
    }

    SECTION("Penalty is symmetric for different gaps") {
        auto p1 = scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_A);
        auto p2 = scorer.shift_penalty(ColumnType::INS_A, ColumnType::DEL_A);
        CHECK(p1 == p2);
    }
}

TEST_CASE("ShiftScoring - Different delta values scale correctly", "[shift_scoring]") {
    SECTION("Delta = 50") {
        ShiftScoring scorer(score_t(50));
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A) == score_t(50));
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_A) == score_t(100));
    }

    SECTION("Delta = 200") {
        ShiftScoring scorer(score_t(200));
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A) == score_t(200));
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_A) == score_t(400));
    }

    SECTION("Delta = 0 (no shift penalty)") {
        ShiftScoring scorer(score_t(0));
        CHECK(scorer.shift_penalty(ColumnType::MATCH, ColumnType::DEL_A) == score_t(0));
        CHECK(scorer.shift_penalty(ColumnType::DEL_A, ColumnType::INS_A) == score_t(0));
    }
}

TEST_CASE("ShiftScoring - Complete 5×5 lookup table verification", "[shift_scoring]") {
    ShiftScoring scorer(score_t(100));

    // Define all column types for iteration
    ColumnType types[] = {
        ColumnType::MATCH,
        ColumnType::DEL_A,
        ColumnType::INS_A,
        ColumnType::DEL_B,
        ColumnType::INS_B
    };

    // Expected penalties for the 5×5 table
    // Rows: c_U, Columns: c_V
    // Order: MATCH, DEL_A, INS_A, DEL_B, INS_B
    int expected[5][5] = {
        //  M    DA   IA   DB   IB
        {   0,  100, 100, 100, 100 },  // MATCH
        { 100,   0,  200,   0, 200 },  // DEL_A
        { 100, 200,   0,  200,   0 },  // INS_A
        { 100,   0,  200,   0, 200 },  // DEL_B
        { 100, 200,   0,  200,   0 }   // INS_B
    };

    SECTION("Verify entire 5×5 table") {
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
                auto penalty = scorer.shift_penalty(types[i], types[j]);
                INFO("c_U = " << i << ", c_V = " << j);
                CHECK(penalty == score_t(expected[i][j]));
            }
        }
    }
}
