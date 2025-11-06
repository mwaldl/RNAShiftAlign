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
#include <LocARNA/stopwatch.hh>

//! Version string (from configure.ac via autoconf system)
const std::string VERSION_STRING = (std::string)PACKAGE_STRING;

int
main(int argc, char **argv) {
    LocARNA::StopWatch stopwatch(false);
    stopwatch.start("total");

    std::cout << VERSION_STRING << std::endl;
    std::cout << "RNAShiftAlign - Sankoff-Style RNA Bi-Alignments" << std::endl;
    std::cout << "Under development..." << std::endl;

    stopwatch.stop("total");
    stopwatch.print_info(std::cerr);

    return 0;
}
