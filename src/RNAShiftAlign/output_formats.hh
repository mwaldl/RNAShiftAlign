/**
 * \file output_formats.hh
 *
 * \brief Structured output for RNAShiftAlign results
 *
 * Defines AlignmentResult and format functions.
 * To add TSV: implement format_result_tsv(const AlignmentResult&).
 */

#ifndef RNASHIFTALIGN_OUTPUT_FORMATS_HH
#define RNASHIFTALIGN_OUTPUT_FORMATS_HH

#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>
#include "RNAShiftAlign/shift_aligner.hh"

namespace RNAShiftAlign {

struct AlignmentResult {
    // Tool metadata
    std::string tool_version;
    std::string timestamp;        ///< ISO 8601, filled at run start

    // Input
    std::string file_a;           ///< Path to input file A, or "<builtin>"
    std::string file_b;           ///< Path to input file B, or "<builtin>"
    std::string seq_a;
    std::string seq_b;

    // Settings used
    ScoringParams scoring_params;
    PreprocessingParams prep_params;

    // Alignment output
    ShiftAligner::score_t score;
    bool shifts_occurred;
    std::string alignment_U_a;        ///< sequence of RNA A; based on layer U
    std::string alignment_U_b;        ///< sequence of RNA B; based on layer U
    std::string alignment_V_a;        ///< sequence of RNA A; based on layer V
    std::string alignment_V_b;        ///< sequence of RNA B; based on layer V
    std::string consensus_sequence;   ///< nucleotide where U_A==U_B, '.' otherwise; based on layer U
    std::string consensus_structure;  ///< dot-bracket of consensus structure; based on V layer
    std::string str_a;                ///< dot-bracket of centroid structure of RNA A (prob>=0.5); based on V
    std::string str_b;                ///< dot-bracket of centroid structure of RNA B (prob>=0.5); based on V
    std::string shift_string;         ///< '=' where U==V column, 'S' where gap patterns differ
};

/// Return current UTC time as ISO 8601 string (e.g. "2026-05-18T14:32:00Z").
inline std::string
current_timestamp() {
    std::time_t t = std::time(nullptr);
    std::tm* utc = std::gmtime(&t);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", utc);
    return buf;
}

/// Escape a string for safe embedding in a JSON value.
inline std::string
json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
        case '"':  out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:   out += c;      break;
        }
    }
    return out;
}

inline std::string
format_result_json(const AlignmentResult& r) {
    std::ostringstream o;
    o << "{\n"
      << "  \"tool\": \"rnashiftalign\",\n"
      << "  \"version\": \"" << json_escape(r.tool_version) << "\",\n"
      << "  \"timestamp\": \"" << json_escape(r.timestamp) << "\",\n"
      << "  \"input\": {\n"
      << "    \"file_a\": \"" << json_escape(r.file_a) << "\",\n"
      << "    \"file_b\": \"" << json_escape(r.file_b) << "\",\n"
      << "    \"seq_a\": \"" << json_escape(r.seq_a) << "\",\n"
      << "    \"seq_b\": \"" << json_escape(r.seq_b) << "\"\n"
      << "  },\n"
      << "  \"settings\": {\n"
      << "    \"match\": " << r.scoring_params.match << ",\n"
      << "    \"mismatch\": " << r.scoring_params.mismatch << ",\n"
      << "    \"indel\": " << r.scoring_params.indel << ",\n"
      << "    \"struct_weight\": " << r.scoring_params.struct_weight << ",\n"
      << "    \"tau\": " << r.scoring_params.tau_factor << ",\n"
      << "    \"delta\": " << r.scoring_params.delta << ",\n"
      << "    \"max_shifts\": " << r.scoring_params.max_shifts << ",\n"
      << "    \"min_prob\": " << r.prep_params.min_prob << ",\n"
      << "    \"max_bps_length_ratio\": " << r.prep_params.max_bps_length_ratio << ",\n"
      << "    \"max_diff_am\": " << r.prep_params.max_diff_am << ",\n"
      << "    \"max_diff_at_am\": " << r.prep_params.max_diff_at_am << "\n"
      << "  },\n"
      << "  \"result\": {\n"
      << "    \"score\": " << r.score << ",\n"
      << "    \"shifts_occurred\": " << (r.shifts_occurred ? "true" : "false") << ",\n"
      << "    \"alignment_U\": {\n"
      << "      \"seq_a\": \"" << json_escape(r.alignment_U_a) << "\",\n"
      << "      \"seq_b\": \"" << json_escape(r.alignment_U_b) << "\"\n"
      << "    },\n"
      << "    \"alignment_V\": {\n"
      << "      \"seq_a\": \"" << json_escape(r.alignment_V_a) << "\",\n"
      << "      \"seq_b\": \"" << json_escape(r.alignment_V_b) << "\"\n"
      << "    },\n"
      << "    \"consensus_sequence\": \"" << json_escape(r.consensus_sequence) << "\",\n"
      << "    \"consensus_structure\": \"" << json_escape(r.consensus_structure) << "\",\n"
      << "    \"str_a\": \"" << json_escape(r.str_a) << "\",\n"
      << "    \"str_b\": \"" << json_escape(r.str_b) << "\",\n"
      << "    \"shift_string\": \"" << json_escape(r.shift_string) << "\"\n"
      << "  }\n"
      << "}\n";
    return o.str();
}

} // namespace RNAShiftAlign

#endif // RNASHIFTALIGN_OUTPUT_FORMATS_HH
