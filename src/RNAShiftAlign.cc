/**
 * \file RNAShiftAlign.cc
 *
 * \brief Main program for RNAShiftAlign
 *
 * RNAShiftAlign: Sankoff-Style RNA Bi-Alignments with Shift Detection
 *
 * Copyright (C) Maria Waldl <maria@bioinf.uni-leipzig.de>
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <LocARNA/stopwatch.hh>
#include "RNAShiftAlign/shift_aligner.hh"

using namespace RNAShiftAlign;

/**
 * @brief Read first sequence from FASTA file
 *
 * Reads a simple FASTA format file and returns the first sequence found.
 * Format: >header\nSEQUENCE (can span multiple lines)
 *
 * @param filename Path to FASTA file
 * @return RNA sequence string (without header or whitespace)
 * @throws std::runtime_error if file cannot be opened or no sequence found
 */
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
            // Header line
            if (in_sequence) {
                // Already read a sequence, return it
                break;
            }
            in_sequence = true;
        } else if (in_sequence) {
            // Sequence line - accumulate
            sequence += line;
        }
    }

    if (sequence.empty()) {
        throw std::runtime_error("No sequence found in file: " + filename);
    }

    return sequence;
}

//! Version string (from configure.ac via autoconf system)
const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

int
main(int argc, char **argv) {
    LocARNA::StopWatch stopwatch(false);
    stopwatch.start("total");

    std::cout << VERSION_STRING << std::endl;
    std::cout << "RNAShiftAlign - Sankoff-Style RNA Bi-Alignments" << std::endl;
    std::cout << "Under development..." << std::endl;
    std::cout << std::endl;

    // ========== Read input sequences ==========
    std::string seqA, seqB;

    if (argc >= 3) {
        // Read from files
        std::string fileA = argv[1];
        std::string fileB = argv[2];

        std::cout << "Reading sequence A from: " << fileA << std::endl;
        std::cout << "Reading sequence B from: " << fileB << std::endl;

        try {
            seqA = read_fasta(fileA);
            seqB = read_fasta(fileB);
        } catch (const std::exception& e) {
            std::cerr << "Error reading input files: " << e.what() << std::endl;
            std::cerr << std::endl;
            std::cerr << "Usage: " << argv[0] << " <seqA.fasta> <seqB.fasta>" << std::endl;
            std::cerr << "  Aligns two RNA sequences from FASTA files" << std::endl;
            std::cerr << std::endl;
            std::cerr << "  Or run without arguments to use default test sequences" << std::endl;
            return 1;
        }
    } else {
        // Use default test sequences
        std::cout << "No input files specified, using default test sequences" << std::endl;
        seqA = "GGGAAACGGGGC";
        seqB = "GGGUAACGGGGCC";
    }

    std::cout << "Sequence A (" << seqA.length() << " nt): " << seqA << std::endl;
    std::cout << "Sequence B (" << seqB.length() << " nt): " << seqB << std::endl;
    std::cout << std::endl;

    // ========== Setup parameters ==========
    PreprocessingParams prep_params;
    prep_params.min_prob = 0.001; // locarna default 0.001
    prep_params.max_diff_am = 30;

    ScoringParams scoring_params;
    scoring_params.match = 50;
    scoring_params.mismatch = 0;
    scoring_params.indel = -150; // open 750 in locarna
    scoring_params.struct_weight = 150;
    scoring_params.delta = -200;
    scoring_params.max_shifts = 2;
    // unpaired-penalty = 0
    // tau = 50

    std::cout << "Parameters:" << std::endl;
    std::cout << "  match=" << scoring_params.match
              << ", mismatch=" << scoring_params.mismatch
              << ", indel=" << scoring_params.indel << std::endl;
    std::cout << "  struct_weight=" << scoring_params.struct_weight
              << ", delta=" << scoring_params.delta
              << ", max_shifts=" << scoring_params.max_shifts << std::endl;
    std::cout << "  min_prob=" << prep_params.min_prob
              << ", max_diff_am=" << prep_params.max_diff_am << std::endl;
    std::cout << std::endl;

    // ========== Run alignment ==========
    std::cout << "Computing alignment..." << std::endl;
    ShiftAligner aligner(seqA, seqB, prep_params, scoring_params);
    auto score = aligner.align();

    std::cout << "Alignment score: " << score << std::endl;
    std::cout << std::endl;

    stopwatch.stop("total");
    stopwatch.print_info(std::cerr);

    return 0;
}
