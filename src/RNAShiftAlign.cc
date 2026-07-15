/**
 * \file RNAShiftAlign.cc
 *
 * \brief Main program for RNAShiftAlign
 *
 * RNAShiftAlign: Sankoff-Style RNA Bi-Alignments with Shift Detection
 *
 * Copyright (C) Maria Waldl <code@waldl.org >
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <getopt.h>
#include <LocARNA/stopwatch.hh>
#include "RNAShiftAlign/shift_aligner.hh"
#include "RNAShiftAlign/output_formats.hh"

using namespace RNAShiftAlign;

std::string
read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    std::string sequence;
    bool in_sequence = false;

    while (std::getline(file, line)) {
        // Remove leading/trailing whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;
        // Skip comment lines
        if (line[0] == '#') continue;
        if (line[0] == '>') {
            // only read first sequence in multi-fasta
            if (in_sequence) break;
            in_sequence = true;
        } else if (in_sequence) {
            sequence += line;
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error("No sequence found in file: " + filename);
    }
    return sequence;
}

const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

static void
print_usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [OPTIONS] <seqA.fasta> <seqB.fasta>\n"
        "       " << prog << " [OPTIONS]          (uses built-in test sequences)\n"
        "\n"
        "Sankoff-style RNA bi-alignment with shift detection.\n"
        "\n"
        "General:\n"
        "  -h, --help                  Show this help and exit\n"
        "  -V, --version               Show version and exit\n"
        "  -v, --verbose LEVEL         Verbosity 0-3 (default: 0)\n"
        "                              0=silent, 1=basic, 2=detailed, 3=debug\n"
        "  -q, --quiet                 Suppress human-readable output on stdout\n"
        "  --output-file FILE          Write result as JSON to FILE\n"
        "\n"
        "Preprocessing:\n"
        "  --min-prob FLOAT            Min base pair probability (default: 0.001)\n"
        "  --max-bps-length-ratio FLOAT  Max base-pairs/length ratio (default: 0.0 = off)\n"
        "  --max-diff-am INT           Max arc length difference (default: 100; -1=off)\n"
        "  --max-diff-at-am INT        Max position difference at arc ends (default: -1=off)\n"
        "\n"
        "Scoring:\n"
        "  --match INT                 Match score (default: 50)\n"
        "  --mismatch INT              Mismatch score (default: 0)\n"
        "  --indel INT                 Gap penalty (default: -100)\n"
        "  --struct-weight INT         Structure weight (default: 150)\n"
        "  --delta INT                 Shift penalty (default: -200)\n"
        "  --max-shifts INT            Max allowed shifts delta_max (default: 1)\n";
}

int
main(int argc, char **argv) {
    LocARNA::StopWatch stopwatch(false);
    stopwatch.start("total");

    // ========== Defaults ==========
    int verbose = 0;
    bool quiet = false;
    std::string output_file;
    std::string timestamp = current_timestamp();

    // Defaults live in the parameter structs (see shift_aligner.hh) and are
    // mirrored in print_usage(). Options below override individual fields only.
    PreprocessingParams prep_params;
    ScoringParams scoring_params;

    // ========== Option table ==========
    enum LongOpts {
        OPT_MIN_PROB = 256,
        OPT_MAX_BPS_LEN_RATIO,
        OPT_MAX_DIFF_AM,
        OPT_MAX_DIFF_AT_AM,
        OPT_MATCH,
        OPT_MISMATCH,
        OPT_INDEL,
        OPT_STRUCT_WEIGHT,
        OPT_DELTA,
        OPT_MAX_SHIFTS,
        OPT_OUTPUT_FILE,
    };

    static struct option long_options[] = {
        {"help",                  no_argument,       nullptr, 'h'},
        {"version",               no_argument,       nullptr, 'V'},
        {"verbose",               required_argument, nullptr, 'v'},
        {"quiet",                 no_argument,       nullptr, 'q'},
        {"output-file",           required_argument, nullptr, OPT_OUTPUT_FILE},
        {"min-prob",              required_argument, nullptr, OPT_MIN_PROB},
        {"max-bps-length-ratio",  required_argument, nullptr, OPT_MAX_BPS_LEN_RATIO},
        {"max-diff-am",           required_argument, nullptr, OPT_MAX_DIFF_AM},
        {"max-diff-at-am",        required_argument, nullptr, OPT_MAX_DIFF_AT_AM},
        {"match",                 required_argument, nullptr, OPT_MATCH},
        {"mismatch",              required_argument, nullptr, OPT_MISMATCH},
        {"indel",                 required_argument, nullptr, OPT_INDEL},
        {"struct-weight",         required_argument, nullptr, OPT_STRUCT_WEIGHT},
        {"delta",                 required_argument, nullptr, OPT_DELTA},
        {"max-shifts",            required_argument, nullptr, OPT_MAX_SHIFTS},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "hVv:q", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            print_usage(argv[0]);
            return 0;
        case 'V':
            std::cout << VERSION_STRING << "\n";
            return 0;
        case 'v':
            verbose = std::stoi(optarg);
            if (verbose < 0 || verbose > 3) {
                std::cerr << "Error: --verbose must be 0-3\n";
                return 1;
            }
            break;
        case 'q':
            quiet = true;
            break;
        case OPT_OUTPUT_FILE:
            output_file = optarg;
            break;
        case OPT_MIN_PROB:
            prep_params.min_prob = std::stod(optarg);
            break;
        case OPT_MAX_BPS_LEN_RATIO:
            prep_params.max_bps_length_ratio = std::stod(optarg);
            break;
        case OPT_MAX_DIFF_AM:
            prep_params.max_diff_am = std::stoi(optarg);
            break;
        case OPT_MAX_DIFF_AT_AM:
            prep_params.max_diff_at_am = std::stoi(optarg);
            break;
        case OPT_MATCH:
            scoring_params.match = std::stoi(optarg);
            break;
        case OPT_MISMATCH:
            scoring_params.mismatch = std::stoi(optarg);
            break;
        case OPT_INDEL:
            scoring_params.indel = std::stoi(optarg);
            break;
        case OPT_STRUCT_WEIGHT:
            scoring_params.struct_weight = std::stoi(optarg);
            break;
        case OPT_DELTA:
            scoring_params.delta = std::stoi(optarg);
            break;
        case OPT_MAX_SHIFTS:
            scoring_params.max_shifts = std::stoi(optarg);
            break;
        default:
            print_usage(argv[0]);
            return 1;
        }
    }

    if (!quiet) {
        std::cout << VERSION_STRING << "\n"
                  << "RNAShiftAlign - Sankoff-Style RNA Bi-Alignments\n"
                  << "Under development...\n\n";
    }

    // ========== Read input sequences ==========
    std::string seqA, seqB;
    std::string fileA_label, fileB_label;
    int remaining = argc - optind;

    if (remaining >= 2) {
        fileA_label = argv[optind];
        fileB_label = argv[optind + 1];
        if (!quiet) {
            std::cout << "Reading sequence A from: " << fileA_label << "\n";
            std::cout << "Reading sequence B from: " << fileB_label << "\n";
        }
        try {
            seqA = read_fasta(fileA_label);
            seqB = read_fasta(fileB_label);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input files: " << e.what() << "\n";
            return 1;
        }
    } else if (remaining == 0) {
        fileA_label = "<builtin>";
        fileB_label = "<builtin>";
        if (!quiet) {
            std::cout << "No input files specified, using default test sequences\n";
        }
        seqA = "CCCCAAAGGG";
        seqB = "CCCAAAUGGG";
    } else {
        std::cerr << "Error: provide either 0 or 2 positional arguments (FASTA files)\n";
        print_usage(argv[0]);
        return 1;
    }

    if (!quiet) {
        std::cout << "Sequence A (" << seqA.length() << " nt): " << seqA << "\n";
        std::cout << "Sequence B (" << seqB.length() << " nt): " << seqB << "\n\n";

        std::cout << "Parameters:\n"
                  << "  match=" << scoring_params.match
                  << ", mismatch=" << scoring_params.mismatch
                  << ", indel=" << scoring_params.indel << "\n"
                  << "  struct_weight=" << scoring_params.struct_weight
                  << ", delta=" << scoring_params.delta
                  << ", max_shifts=" << scoring_params.max_shifts << "\n"
                  << "  min_prob=" << prep_params.min_prob
                  << ", max_bps_length_ratio=" << prep_params.max_bps_length_ratio
                  << ", max_diff_am=" << prep_params.max_diff_am
                  << ", max_diff_at_am=" << prep_params.max_diff_at_am << "\n\n";
    }

    // ========== Run alignment ==========
    if (!quiet) std::cout << "Computing alignment...\n";
    ShiftAligner aligner(seqA, seqB, prep_params, scoring_params, verbose);
    auto score = aligner.align();

    if (!quiet) std::cout << "Alignment score: " << score << "\n\n";

    if (!quiet) std::cout << "Performing traceback...\n";
    aligner.traceback();

    if (!quiet) {
        std::cout << "\n";
        std::cout << aligner.format_alignment() << "\n";
    }

    // ========== JSON output ==========
    if (!output_file.empty()) {
        auto [U_a, U_b] = aligner.get_alignment_U();
        auto [V_a, V_b] = aligner.get_alignment_V();

        bool shifts_occurred = false;
        if (U_a.size() == V_a.size()) {
            for (size_t i = 0; i < U_a.size(); ++i) {
                if ((U_a[i] == '-') != (V_a[i] == '-') ||
                    (U_b[i] == '-') != (V_b[i] == '-')) {
                    shifts_occurred = true;
                    break;
                }
            }
        }

        AlignmentResult result;
        result.tool_version  = VERSION_STRING;
        result.timestamp     = timestamp;
        result.file_a        = fileA_label;
        result.file_b        = fileB_label;
        result.seq_a         = seqA;
        result.seq_b         = seqB;
        result.scoring_params = scoring_params;
        result.prep_params   = prep_params;
        result.score              = score;
        result.shifts_occurred    = shifts_occurred;
        result.alignment_U_a      = U_a;
        result.alignment_U_b      = U_b;
        result.alignment_V_a      = V_a;
        result.alignment_V_b      = V_b;
        result.consensus_sequence  = aligner.get_consensus_sequence();
        result.consensus_structure = aligner.get_consensus_structure();
        result.str_a               = aligner.get_dot_bracket_A();
        result.str_b               = aligner.get_dot_bracket_B();
        result.shift_string        = aligner.get_shift_string();

        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: cannot open output file: " << output_file << "\n";
            return 1;
        }
        out << format_result_json(result);
        if (!quiet) std::cout << "JSON result written to: " << output_file << "\n";
    }

    stopwatch.stop("total");
    stopwatch.print_info(std::cerr);

    return 0;
}
