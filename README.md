# RNAShiftAlign

**Sankoff-Style RNA Bi-Alignments with Shift Detection**

RNAShiftAlign simultaneously aligns and folds two RNA sequences, and detects
*shifts* — cases where the conserved structure has moved locally relative to
the homologous sequence positions over the course of evolution.

## Overview

RNA molecules can evolve under independent selection on their sequence and on
their secondary structure. This produces cases where:

- structurally equivalent base pairs involve non-homologous sequence positions, and
- homologous sequence positions no longer preserve a base pair.

A standard alignment must choose one column assignment; it cannot express that
the sequence-optimal and structure-optimal alignments disagree. RNAShiftAlign
computes a *bi-alignment*: two coupled alignments of the same sequence pair —
a **sequence layer** (U) and a **structure layer** (V) — together with a
per-column **shift annotation** (W) that marks where the two layers diverge.
The number of positions the two layers may drift apart is bounded by
`--max-shifts`.

The method is an extension of the Sankoff algorithm for simultaneous alignment
and folding, described in *"Incongruences Between Sequence and Secondary
Structure Alignments of Nucleic Acids"* (Waldl et al. 2026).

## Dependencies

- A **C++14** compiler
- **[LocARNA](https://github.com/s-will/LocARNA) ≥ 2.0** 
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) ≥ 2.5.1** — pulled in
  automatically as a dependency of LocARNA
- **pkg-config** and the **GNU Autotools** (autoconf, automake)

LocARNA's `pkg-config` file (`LocARNA-2.0.pc`) must be discoverable, i.e. its
directory must be on `PKG_CONFIG_PATH`. LocARNA and ViennaRNA can be installed,
for example, via conda/bioconda:

```bash
conda install -c bioconda locarna viennarna
```

## Installation

Standard Autotools, out-of-tree build:

```bash
autoreconf -i                       # generate ./configure (first checkout only)
./configure                         # or: ./configure --prefix=$HOME/.local
make
make install                        # optional; installs the rnashiftalign binary
```

If `configure` cannot find LocARNA, point `pkg-config` at it:

```bash
PKG_CONFIG_PATH=/path/to/locarna/lib/pkgconfig ./configure
```

## Usage

```
rnashiftalign [OPTIONS] <seqA.fasta> <seqB.fasta>
rnashiftalign [OPTIONS]          # no files → uses a small built-in example
```

Each input is a FASTA file; only the first sequence in a multi-FASTA file is
read. With no positional arguments, the tool runs on a built-in example pair so
you can check your build quickly:

```console
$ rnashiftalign
...
4-Way Bi-Alignment:
U_A:  CCCCAAAGGG
Cseq: CCC.AA.GGG
U_B:  CCCAAAUGGG
Shft: ==========
V_A:  CCCCAAAGGG
db_A: (((....)))
Cons: (((....)))
db_B: (((....)))
V_B:  CCCAAAUGGG
```

Reading the 4-way output:

| Row  | Meaning |
|------|---------|
| `U_A` / `U_B` | the two sequences aligned on the **sequence** layer (U) |
| `Cseq` | consensus sequence (nucleotide where U columns match, `.` otherwise) |
| `Shft` | shift track: `=` where U and V agree, `S` where they shift |
| `V_A` / `V_B` | the two sequences aligned on the **structure** layer (V) |
| `db_A` / `db_B` | each sequence's structure (centroid, base pairs with prob ≥ 0.5), aligned to V |
| `Cons` | consensus structure over the structure layer |

When there are no shifts, the U and V layers are identical and `Shft` is all
`=`. Shifts appear as `S` in the shift track, in columns where the `U_*` and
`V_*` alignments differ. The consensus structure is derived from the base-pair
probabilities in each sequence's structure ensemble, as computed by
RNAfold/ViennaRNA.

### Machine-readable output

`--output-file result.json` writes the full result — inputs, settings,
both alignment layers, consensus sequence/structure, per-sequence structures,
and the shift string — as JSON. Combine with `-q` to suppress the
human-readable stdout report.

### Options

```
General:
  -h, --help                    Show help and exit
  -V, --version                 Show version and exit
  -v, --verbose LEVEL           Verbosity 0-3 (default: 0)
  -q, --quiet                   Suppress human-readable stdout
  --output-file FILE            Write result as JSON to FILE

Preprocessing (base-pair probability filtering):
  --min-prob FLOAT              Min base-pair probability (default: 0.001)
  --max-bps-length-ratio FLOAT  Max base-pairs/length ratio (default: 0.0 = off)
  --max-diff-am INT             Max arc-length difference (default: 100; -1 = off)
  --max-diff-at-am INT          Max position difference at arc ends (default: -1 = off)

Scoring:
  --match INT                   Match score (default: 50)
  --mismatch INT                Mismatch score (default: 0)
  --indel INT                   Gap penalty (default: -100)
  --struct-weight INT           Structure weight (default: 150)
  --delta INT                   Shift penalty (default: -200)
  --max-shifts INT              Max allowed shifts δ_max (default: 1)
```

The two parameters most worth tuning are `--max-shifts` (how far the sequence
and structure alignments may diverge; `0` forces a single classical alignment)
and `--delta` (how strongly each shift is penalised).

The computation can be sped up by pruning the set of base pairs and arc matches
considered: raise `--min-prob` to require a higher ensemble probability before a
base pair is kept, and use `--max-bps-length-ratio`, `--max-diff-am`, and
`--max-diff-at-am` to restrict which base pairs may be matched.

## How it works (program flow)

1. **Read input** — parse the two FASTA sequences (`RNAShiftAlign.cc`).
2. **Preprocess with LocARNA** — fold each sequence into a base-pair
   probability matrix, keep pairs above `--min-prob`, and enumerate the
   candidate *arc matches* (pairs of base pairs, one from each sequence) that
   the structure alignment may use. The `--max-diff-*` and
   `--max-bps-length-ratio` options prune this candidate set.
3. **Forward pass** — dynamic programming over two coupled matrices
   (`shift_aligner.cc`):
   - the 4-D matrix **M** aligns the sequence layer (U) and structure layer (V)
     jointly, with the constraint that the two layers stay within `--max-shifts`
     columns of each other;
   - the 6-D matrix **D** scores the interior of each arc match, so that a
     matched base pair contributes the alignment of everything nested inside it.

   The two are mutually recursive: M consults D at arc matches, D fills its
   interior via local instances of M. Each column contributes three additive terms:
   - **sequence** — over the U layer: `--match`/`--mismatch` for an aligned pair
     of nucleotides, `--indel` for a gap;
   - **structure** — over the V layer: 0 for an unpaired column. Base pairs are
     not scored column-by-column; instead a whole arc match is rewarded once, by
     an amount derived from the two base pairs' ensemble probabilities and scaled
     by `--struct-weight`;
   - **shift** — `--delta` for each sequence whose gap pattern differs between
     the U and V layers in this column (so 0, one, or two times `--delta`).
4. **Traceback** — recover the optimal bi-alignment from the (re-)filled matrices,
   producing the U and V alignments and the shift annotation.
5. **Output** — assemble the consensus sequence and structure and print the
   4-way alignment; optionally write the JSON result.

## Authors

Maria Waldl and collaborators.

The implementation builds on the [LocARNA](https://github.com/s-will/LocARNA)
library by Sebastian Will.
