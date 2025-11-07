#include "shift_aligner.hh"

#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/pfold_params.hh>

namespace RNAShiftAlign {

ShiftAligner::ShiftAligner(const std::string &seqA,
                           const std::string &seqB,
                           score_t delta,
                           size_type max_shifts)
    : shift_scoring_(delta),
      max_shifts_(max_shifts),
      alignment_score_(score_t(0))
{
    // Store sequences
    seqA_ = std::make_unique<LocARNA::Sequence>("seqA", seqA);
    seqB_ = std::make_unique<LocARNA::Sequence>("seqB", seqB);

    // TODO: Initialize RnaData and Scoring when implementing actual recursion.
    // RnaData construction requires proper LocARNA setup (partition function, etc.)
    // Scoring requires ArcMatches and MatchProbs for structure matching.
    // For skeleton testing, we only need the sequences and DP matrix.

    // Allocate main DP matrix M(y1, y2, y3, y4)
    // Dimensions: (len_A + 1) × (len_B + 1) × shifts × shifts
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    M_ = std::make_unique<LocARNA::ShiftMatrixM<score_t>>(
        lenA + 1,  // Include position 0 for empty prefixes
        lenB + 1,
        max_shifts
    );

    // Initialize M matrix with -infinity (maximization problem)
    M_->fill(LocARNA::infty_score_t::neg_infty);
    // Base case: M(0, 0, 0, 0) = 0 (empty alignment)
    M_->set(0, 0, 0, 0, score_t(0));
}

ShiftAligner::score_t
ShiftAligner::align() {
    // Fill M matrix for unpaired-only alignment
    fill_M_unpaired();

    // Optimal score is at M(len_A, len_B, len_A, len_B)
    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    alignment_score_ = M_->get(lenA, lenB, lenA, lenB);
    return alignment_score_;
}

void
ShiftAligner::traceback() {
    // TODO: Implement traceback in future phase
    // Will reconstruct U, V, W alignments from M matrix
}

void
ShiftAligner::fill_M_unpaired() {
    // TODO: Implement unpaired forward recursion in next phase
    // This will fill M matrix according to Case 1 of Equation 17
}

} // namespace RNAShiftAlign
