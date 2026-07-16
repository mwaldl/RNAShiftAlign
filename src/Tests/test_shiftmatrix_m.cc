#include "catch.hpp"

#include <cassert>
#include <algorithm>
#include <ShiftAlign/shiftmatrix_m.hh>

using namespace RNAShiftAlign;

/**
 * @file test_shiftmatrix_m.cc
 * @brief Unit tests for ShiftMatrixM class
 *
 * Tests the 4D matrix used for bi-alignment dynamic programming,
 * including index validation, maxshift constraint enforcement,
 * and basic operations (fill, transform, get/set).
 */

/**
 * @brief Helper functor for testing transform() method
 */
struct mul2 {
    size_t
    operator()(size_t x) {
        return 2 * x;
    }
};

TEST_CASE("ShiftMatrixM - Basic Operations") {
    size_t a = 3;
    size_t b = 4;
    size_t max_shift = 1;

    ShiftMatrixM<size_t> m(a, b, max_shift);

    SECTION("Can be resized, filled, and transformed") {
        m.resize(a, b + 1, max_shift);

        m.fill(2);
        m.set(1, 2, 1, 3, 5);

        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < b; j++) {
                m(i, j, i, j) += i + 2 * j;
            }
        }
        m.transform(mul2());

        REQUIRE(m.get(1, 2, 1, 3) == 10);

        bool reread_ok = true;
        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < b; j++) {
                reread_ok &= (m(i, j, i, j) == 2 * (2 + i + 2 * j));
            }
        }
        REQUIRE(reread_ok);
    }

    SECTION("Accepts valid indices with no shift") {
        CHECK_NOTHROW(m.get(0, 0, 0, 0));
        CHECK_NOTHROW(m.get(1, 1, 1, 1));
        CHECK_NOTHROW(m.get(2, 3, 2, 3));
    }

    SECTION("Accepts valid indices with maxshift") {
        CHECK_NOTHROW(m.get(1, 3, 2, 3));  // shift_a = 1, shift_b = 0
        CHECK_NOTHROW(m.get(0, 3, 1, 3));  // shift_a = 1, shift_b = 0
        CHECK_NOTHROW(m.get(1, 1, 1, 2));  // shift_a = 0, shift_b = 1
        CHECK_NOTHROW(m.get(2, 2, 1, 1));  // shift_a = 1, shift_b = 1
    }
}

TEST_CASE("ShiftMatrixM - Index Validation") {
    size_t a = 3;
    size_t b = 4;
    size_t max_shift = 1;

    ShiftMatrixM<size_t> m(a, b, max_shift);

    SECTION("Throws exception when y1 index is out of bounds") {
        CHECK_THROWS(m.get(3, 0, 3, 0));   // y1 = adim
        CHECK_THROWS(m.get(100, 0, 100, 0));
    }

    SECTION("Throws exception when y2 index is out of bounds") {
        CHECK_THROWS(m.get(0, 4, 0, 4));   // y2 = bdim
        CHECK_THROWS(m.get(0, 100, 0, 100));
    }

    SECTION("Throws exception when y3 index is out of bounds") {
        CHECK_THROWS(m.get(0, 0, 3, 0));   // y3 = adim
        CHECK_THROWS(m.get(0, 0, 100, 0));
    }

    SECTION("Throws exception when y4 index is out of bounds") {
        CHECK_THROWS(m.get(0, 0, 0, 4));   // y4 = bdim
        CHECK_THROWS(m.get(0, 0, 0, 100));
    }
}

TEST_CASE("ShiftMatrixM - Maxshift Constraint") {
    size_t a = 5;
    size_t b = 5;
    size_t max_shift = 1;

    ShiftMatrixM<size_t> m(a, b, max_shift);

    SECTION("Rejects shift_a > maxshift") {
        CHECK_THROWS(m.get(0, 0, 2, 0));  // |y1-y3| = 2 > maxshift
        CHECK_THROWS(m.get(2, 0, 0, 0));  // |y1-y3| = 2 > maxshift
    }

    SECTION("Rejects shift_b > maxshift") {
        CHECK_THROWS(m.get(0, 0, 0, 2));  // |y2-y4| = 2 > maxshift
        CHECK_THROWS(m.get(0, 2, 0, 0));  // |y2-y4| = 2 > maxshift
    }

    SECTION("Accepts shift exactly at maxshift") {
        CHECK_NOTHROW(m.get(1, 1, 2, 2));  // shift_a = 1, shift_b = 1
        CHECK_NOTHROW(m.get(2, 2, 1, 1));  // shift_a = 1, shift_b = 1
    }
}

TEST_CASE("ShiftMatrixM - Different maxshift values") {
    SECTION("maxshift = 0 (no shifts allowed)") {
        ShiftMatrixM<int> m(3, 3, 0);

        CHECK_NOTHROW(m.get(0, 0, 0, 0));
        CHECK_NOTHROW(m.get(1, 1, 1, 1));
        CHECK_THROWS(m.get(0, 0, 1, 0));  // any shift forbidden
        CHECK_THROWS(m.get(0, 0, 0, 1));
    }

    SECTION("maxshift = 2 (larger shifts allowed)") {
        ShiftMatrixM<int> m(5, 5, 2);

        CHECK_NOTHROW(m.get(0, 0, 2, 2));  // shift = 2
        CHECK_NOTHROW(m.get(2, 2, 0, 0));  // shift = 2
        CHECK_THROWS(m.get(0, 0, 3, 0));   // shift = 3 > maxshift
    }
}

TEST_CASE("ShiftMatrixM - Clear operation") {
    ShiftMatrixM<int> m(3, 4, 1);
    m.fill(42);

    m.clear();

    auto sizes = m.sizes();
    REQUIRE(std::get<0>(sizes) == 0);  // adim
    REQUIRE(std::get<1>(sizes) == 0);  // bdim
    REQUIRE(m.get_maxshift() == 0);
}

TEST_CASE("ShiftMatrixM - Empty constructor") {
    ShiftMatrixM<double> m;

    REQUIRE(m.get_maxshift() == 0);
    auto sizes = m.sizes();
    REQUIRE(std::get<0>(sizes) == 0);
    REQUIRE(std::get<1>(sizes) == 0);
}
