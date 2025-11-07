#include "shift_aligner.hh"

#include <LocARNA/rna_ensemble.hh>
#include <LocARNA/pfold_params.hh>
#include <LocARNA/anchor_constraints.hh>
#include <LocARNA/trace_controller.hh>
#include <LocARNA/multiple_alignment.hh>

namespace RNAShiftAlign {

ShiftAligner::ShiftAligner(const std::string &seqA,
                           const std::string &seqB,
                           const PreprocessingParams &prep_params,
                           const ScoringParams &scoring_params)
    : shift_scoring_(score_t(scoring_params.delta)),
      max_shifts_(scoring_params.max_shifts),
      alignment_score_(score_t(0))
{
    // Step 1: Create sequences
    seqA_ = std::make_unique<LocARNA::Sequence>("seqA", seqA);
    seqB_ = std::make_unique<LocARNA::Sequence>("seqB", seqB);

    size_type lenA = seqA_->length();
    size_type lenB = seqB_->length();

    // Step 2: Compute base pair probabilities via ViennaRNA partition function
    // Sequence → MultipleAlignment → RnaEnsemble → RnaData
    LocARNA::PFoldParams pfold_params;  // Use default ViennaRNA parameters

    // Create MultipleAlignment from single sequences
    auto ma_A = std::make_unique<LocARNA::MultipleAlignment>("seqA", seqA);
    auto ma_B = std::make_unique<LocARNA::MultipleAlignment>("seqB", seqB);

    // Create RnaEnsemble from MultipleAlignments
    auto rna_ensembleA = std::make_unique<LocARNA::RnaEnsemble>(
        *ma_A,
        pfold_params,
        false,  // no stacking
        false   // not in_loop
    );

    auto rna_ensembleB = std::make_unique<LocARNA::RnaEnsemble>(
        *ma_B,
        pfold_params,
        false,  // no stacking
        false   // not in_loop
    );

    // Create RnaData from ensembles
    rna_dataA_ = std::make_unique<LocARNA::RnaData>(
        *rna_ensembleA,
        prep_params.min_prob,
        prep_params.max_bps_length_ratio,
        pfold_params
    );

    rna_dataB_ = std::make_unique<LocARNA::RnaData>(
        *rna_ensembleB,
        prep_params.min_prob,
        prep_params.max_bps_length_ratio,
        pfold_params
    );

    // Step 3: Compute arc matches between the two sequences
    // Need to provide TraceController and AnchorConstraints (use defaults)

    // Default constraints: no anchors
    LocARNA::AnchorConstraints seq_constraints(lenA, "", lenB, "", true);

    // Default trace controller: no restrictions (use -1 for max_diff to disable, false for no relaxation)
    auto trace_controller = std::make_unique<LocARNA::TraceController>(
        *seqA_,
        *seqB_,
        nullptr,  // no reference alignment
        -1,       // max_diff = -1 (no restriction)
        false     // no relaxation
    );

    // Compute arc matches with proper parameters
    size_type max_diff_am_value = (prep_params.max_diff_am != -1)
        ? size_type(prep_params.max_diff_am)
        : std::max(lenA, lenB);
    size_type max_diff_at_am_value = (prep_params.max_diff_at_am != -1)
        ? size_type(prep_params.max_diff_at_am)
        : std::max(lenA, lenB);

    arc_matches_ = std::make_unique<LocARNA::ArcMatches>(
        *rna_dataA_,
        *rna_dataB_,
        prep_params.min_prob,
        max_diff_am_value,
        max_diff_at_am_value,
        *trace_controller,
        seq_constraints
    );

    // Step 4: Initialize scoring scheme using named argument pattern
    // Compute expected probabilities
    double exp_probA = rna_dataA_->arc_cutoff_prob();
    double exp_probB = rna_dataB_->arc_cutoff_prob();

    // Create scoring parameters using named argument pattern
    auto locarna_scoring_params = LocARNA::ScoringParams(
        LocARNA::ScoringParams::match(scoring_params.match),
        LocARNA::ScoringParams::mismatch(scoring_params.mismatch),
        LocARNA::ScoringParams::indel(scoring_params.indel),
        LocARNA::ScoringParams::indel_opening(scoring_params.indel_opening),
        LocARNA::ScoringParams::struct_weight(scoring_params.struct_weight),
        LocARNA::ScoringParams::tau_factor(scoring_params.tau_factor),
        LocARNA::ScoringParams::exp_probA(exp_probA),
        LocARNA::ScoringParams::exp_probB(exp_probB)
    );

    // Initialize LocARNA Scoring object (pass nullptr for match_probs - not needed for basic scoring)
    locarna_scoring_ = std::make_unique<LocARNA::Scoring>(
        *seqA_,
        *seqB_,
        *rna_dataA_,
        *rna_dataB_,
        *arc_matches_,
        nullptr,  // match_probs not needed for basic scoring
        locarna_scoring_params
    );

    // Step 5: Allocate shift-aware DP matrix M(y1, y2, y3, y4)
    M_ = std::make_unique<LocARNA::ShiftMatrixM<score_t>>(
        lenA + 1,
        lenB + 1,
        max_shifts_
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
