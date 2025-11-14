#include "catch.hpp"

#include <cassert>
#include <algorithm>
#include <LocARNA/scoring_fwd.hh>
#include <RNAShiftAlign/shiftmatrix_d.hh>

using namespace LocARNA;

/** @file Unit tests for ShiftOffsetMatrix and ShiftMatrixD classes
 *
 * These classes implement the 6D D matrix for shift-aware structure alignment.
 * ShiftOffsetMatrix: 4D matrix storing shifts at arc endpoints
 * ShiftMatrixD: 6D matrix = 2D base pair matrix × 4D offset matrix
 */

TEST_CASE("ShiftOffsetMatrix - Basic operations", "[shiftoffsetmatrix]") {
    size_t max_shift = 1;
    ShiftOffsetMatrix<score_t> offset_matrix(max_shift);

    SECTION("Construction and sizing") {
        REQUIRE(offset_matrix.get_maxshift() == max_shift);

        auto sizes = offset_matrix.sizes();
        size_t expected_dim = 2 * max_shift + 1;  // 3 for maxshift=1
        REQUIRE(std::get<0>(sizes) == expected_dim);
        REQUIRE(std::get<1>(sizes) == expected_dim);
        REQUIRE(std::get<2>(sizes) == expected_dim);
        REQUIRE(std::get<3>(sizes) == expected_dim);
    }

    SECTION("Set and get with valid shifts") {
        // Test with shifts in range [-1, 1]
        score_t test_value = 42;
        offset_matrix.set(0, 0, 0, 0, test_value);
        REQUIRE(offset_matrix.get(0, 0, 0, 0) == test_value);

        offset_matrix.set(-1, 1, 0, -1, 100);
        REQUIRE(offset_matrix.get(-1, 1, 0, -1) == 100);

        offset_matrix.set(1, 1, 1, 1, 200);
        REQUIRE(offset_matrix.get(1, 1, 1, 1) == 200);
    }

    SECTION("Operator() access") {
        offset_matrix(0, 0, 0, 0) = 55;
        REQUIRE(offset_matrix(0, 0, 0, 0) == 55);

        offset_matrix(-1, -1, 1, 1) = 77;
        REQUIRE(offset_matrix(-1, -1, 1, 1) == 77);
    }

    SECTION("Fill operation") {
        offset_matrix.fill(99);
        // Check a few positions
        REQUIRE(offset_matrix(0, 0, 0, 0) == 99);
        REQUIRE(offset_matrix(-1, 0, 1, 0) == 99);
        REQUIRE(offset_matrix(1, 1, -1, -1) == 99);
    }

    SECTION("Resize operation") {
        offset_matrix.resize(2);  // Increase maxshift to 2
        REQUIRE(offset_matrix.get_maxshift() == 2);

        auto sizes = offset_matrix.sizes();
        size_t expected_dim = 5;  // 2*2+1
        REQUIRE(std::get<0>(sizes) == expected_dim);

        // Should be able to access with larger shifts
        CHECK_NOTHROW(offset_matrix.set(-2, 2, -2, 2, 123));
        REQUIRE(offset_matrix.get(-2, 2, -2, 2) == 123);
    }
}

TEST_CASE("ShiftMatrixD - 6D matrix operations", "[shiftmatrixd]") {
    size_t a_bps = 3;  // Number of base pairs in RNA A
    size_t b_bps = 4;  // Number of base pairs in RNA B
    size_t max_shift = 1;

    typedef BasePairs__Arc Arc;

    // Create arcs: Arc(index, left_pos, right_pos)
    Arc arc_a(1, 1, 5);   // Base pair at positions 1-5 in RNA A
    Arc arc_b(3, 2, 7);   // Base pair at positions 2-7 in RNA B

    ShiftMatrixD<score_t> d(a_bps, b_bps, max_shift);

    SECTION("Construction and sizing") {
        REQUIRE(d.get_maxshift() == max_shift);

        auto sizes = d.sizes();
        REQUIRE(std::get<0>(sizes) == a_bps);
        REQUIRE(std::get<1>(sizes) == b_bps);
    }

    SECTION("Create offset matrix and access elements") {
        // Must create offset matrix before accessing
        d.create_offsetmatrix(1, 3);  // For base pair indices (1,3)

        // Set a value
        score_t test_val = 10;
        d.set(arc_a, arc_b, 1, 2, 5, 7, test_val);

        // Get should return the same value
        REQUIRE(d.get(arc_a, arc_b, 1, 2, 5, 7) == test_val);
    }

    SECTION("Multiple offset matrices") {
        // Create offset matrices for different bp pairs
        d.create_offsetmatrix(0, 0);
        d.create_offsetmatrix(1, 3);
        d.create_offsetmatrix(2, 1);

        // Set different values in each
        d.set(Arc(0, 0, 3), Arc(0, 0, 4), 0, 0, 3, 4, 100);
        d.set(arc_a, arc_b, 1, 2, 5, 7, 200);
        d.set(Arc(2, 5, 10), Arc(1, 3, 8), 5, 3, 10, 8, 300);

        // Verify independent storage
        REQUIRE(d.get(Arc(0, 0, 3), Arc(0, 0, 4), 0, 0, 3, 4) == 100);
        REQUIRE(d.get(arc_a, arc_b, 1, 2, 5, 7) == 200);
        REQUIRE(d.get(Arc(2, 5, 10), Arc(1, 3, 8), 5, 3, 10, 8) == 300);
    }

    SECTION("Boundary shifts") {
        d.create_offsetmatrix(1, 3);

        // Test at shift boundaries (-1, 0, 1)
        // arc_b.right() = 7, so valid y2 values are 6, 7, 8 (shifts -1, 0, +1)
        CHECK_NOTHROW(d.set(arc_a, arc_b, 1, 2, 5, 6, 50));  // y2 at boundary (-1)
        CHECK_NOTHROW(d.set(arc_a, arc_b, 1, 2, 5, 7, 60));  // y2 at boundary (0)
        CHECK_NOTHROW(d.set(arc_a, arc_b, 1, 2, 5, 8, 70));  // y2 at boundary (+1)

        REQUIRE(d.get(arc_a, arc_b, 1, 2, 5, 6) == 50);
        REQUIRE(d.get(arc_a, arc_b, 1, 2, 5, 7) == 60);
        REQUIRE(d.get(arc_a, arc_b, 1, 2, 5, 8) == 70);
    }

    SECTION("Accepts various valid arc configurations") {
        d.create_offsetmatrix(1, 3);

        // Different combinations within shift range
        // arc_a.left()=1, arc_a.right()=5, arc_b.left()=2, arc_b.right()=7
        // Valid x1: 0,1,2 (shifts -1,0,+1 from arc_a.left()=1)
        // Valid x2: 1,2,3 (shifts -1,0,+1 from arc_b.left()=2)
        // Valid y1: 4,5,6 (shifts -1,0,+1 from arc_a.right()=5)
        // Valid y2: 6,7,8 (shifts -1,0,+1 from arc_b.right()=7)
        CHECK_NOTHROW(d.get(arc_a, arc_b, 1, 2, 5, 7));  // All zero shifts
        CHECK_NOTHROW(d.get(arc_a, arc_b, 0, 1, 4, 6));  // All -1 shifts
        CHECK_NOTHROW(d.get(arc_a, arc_b, 2, 3, 6, 8));  // All +1 shifts
        CHECK_NOTHROW(d.get(arc_a, arc_b, 1, 3, 4, 8));  // Mixed shifts
    }
}

TEST_CASE("ShiftOffsetMatrix - Addressing correctness", "[shiftoffsetmatrix][addressing]") {
    // Test that addressing formula is correct for 4D matrix
    size_t max_shift = 1;
    ShiftOffsetMatrix<int> matrix(max_shift);

    SECTION("Unique addresses for all valid shift combinations") {
        // For maxshift=1, we have 3^4 = 81 unique positions
        // Fill with unique values and verify no collisions
        int counter = 0;
        for (int x1 = -1; x1 <= 1; ++x1) {
            for (int x2 = -1; x2 <= 1; ++x2) {
                for (int y1 = -1; y1 <= 1; ++y1) {
                    for (int y2 = -1; y2 <= 1; ++y2) {
                        matrix(x1, x2, y1, y2) = counter++;
                    }
                }
            }
        }

        // Verify all values are distinct
        REQUIRE(counter == 81);  // 3^4

        // Check a few specific positions
        REQUIRE(matrix(-1, -1, -1, -1) == 0);
        REQUIRE(matrix(1, 1, 1, 1) == 80);
        REQUIRE(matrix(0, 0, 0, 0) == 40);  // Middle position
    }
}
