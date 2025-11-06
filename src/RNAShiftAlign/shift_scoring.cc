#include "shift_scoring.hh"

namespace RNAShiftAlign {

ShiftScoring::ShiftScoring(score_t delta)
    : delta_(delta) {
}

ShiftScoring::score_t
ShiftScoring::shift_penalty(ColumnType c_U, ColumnType c_V) const {
    // Case 1: Identical column types → no penalty
    // This covers the diagonal of the 15×15 table
    if (c_U == c_V) {
        return score_t(0);
    }

    bool u_is_gap = is_gap(c_U);
    bool v_is_gap = is_gap(c_V);

    // Case 2: One MATCH, one gap → penalty Δ
    // This covers the cross-sections: MATCH vs {DEL_A, INS_A, DEL_B, INS_B}
    if (u_is_gap != v_is_gap) {
        return delta_;
    }

    // Case 3: Both are gaps
    if (u_is_gap && v_is_gap) {
        // Same gap type (both DEL or both INS) → no penalty
        if (gap_types_match(c_U, c_V)) {
            return score_t(0);
        }
        // Different gap types (one DEL, one INS) → penalty 2Δ
        return delta_ + delta_;
    }

    // Should never reach here
    return score_t(0);
}

bool
ShiftScoring::is_gap(ColumnType c) {
    return c != ColumnType::MATCH;
}

bool
ShiftScoring::gap_types_match(ColumnType c_U, ColumnType c_V) {
    // Both are DEL (gap in second sequence)
    bool both_del = (c_U == ColumnType::DEL_A || c_U == ColumnType::DEL_B) &&
                    (c_V == ColumnType::DEL_A || c_V == ColumnType::DEL_B);

    // Both are INS (gap in first sequence)
    bool both_ins = (c_U == ColumnType::INS_A || c_U == ColumnType::INS_B) &&
                    (c_V == ColumnType::INS_A || c_V == ColumnType::INS_B);

    return both_del || both_ins;
}

} // namespace RNAShiftAlign
