# RNAShiftAlign

**Sankoff-Style RNA Bi-Alignments with Shift Detection**

RNAShiftAlign simultaneously aligns and folds two RNA sequences, and detects
*shifts* — cases where the conserved structure has moved locally relative to
the homologous sequence positions over the course of evolution.

## Overview

RNA molecules can evolve under independent selection on their sequence and on
their secondary structure. This can produces cases where:

- structurally equivalent base pairs involve non-homologous sequence positions, and
- homologous sequence positions no longer preserve a base pair.

A standard 2 way alignment must choose one assignment for seqeunce and structure
; it cannot express cases where
the sequence-optimal and structure-optimal alignments disagree. RNAShiftAlign
computes a *bi-alignment*: two coupled alignments of the same RNA pair —
a **sequence layer** (U) and a **structure layer** (V) — together with a
per-column **shift annotation** (W) that marks where the two layers diverge.
The number of positions the two layers may drift apart is bounded by
`--max-shifts`.

The method is an extension of the Sankoff algorithm for simultaneous alignment
and folding, and it is described in *"Incongruences Between Sequence and Secondary
Structure Alignments of Nucleic Acids"* (Waldl et al. 2026).

## Dependencies

- A **C++14** compiler (GCC ≥ 5 or Clang ≥ 3.4)
- **[LocARNA](https://github.com/s-will/LocARNA) ≥ 2.0** — including its
  development files (headers and `LocARNA-2.0.pc`)
- **[ViennaRNA](https://www.tbi.univie.ac.at/RNA/) ≥ 2.5.1** — pulled in
  automatically as a dependency of LocARNA
- **pkg-config** and the **GNU Autotools** (autoconf, automake)

## Installation

### With conda

```bash
conda create -n shift
conda activate shift
conda install -c conda-forge -c bioconda cxx-compiler pkg-config locarna
```

```bash
autoreconf -i                     # generate ./configure (first checkout only)
./configure --prefix=$CONDA_PREFIX
make
make install
```

### Without conda

Install a C++14 compiler, pkg-config and the Autotools, for example, via your system package
manager.

Ubuntu:
```bash
sudo apt install g++ pkg-config autoconf automake
```

LocARNA (and ViennaRNA) must then be built and installed from source; see the
[LocARNA documentation](https://github.com/s-will/LocARNA). If they are
installed under a non-standard prefix, point pkg-config at them:

```bash
export PKG_CONFIG_PATH=/path/to/locarna/lib/pkgconfig
```

```bash
autoreconf -i
./configure --prefix=$HOME/.local
make
make install
```

### Install prefix

`./configure` defaults to `--prefix=/usr/local`, so a plain `make install`
requires root. To install without `sudo`, pass a prefix you own — for example
`$CONDA_PREFIX` inside an active conda environment, or `$HOME/.local` (make sure
`$HOME/.local/bin` is on your `PATH`). Installation can also be skipped entirely
by running the binary from the build tree:

```bash
./src/rnashiftalign
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
U_A:  -AAGGCUCUAUUAACUGGUAUCGGCUAUAG
Cseq: .AA.G.UCUAU.AACUG.UAUC.G.U.UAG
U_B:  -AAUGAUCUAUGAACUGUUAUCUGAUUUAG
Shft: S=========================S===
V_A:  -AAGGCUCUAUUAACUGGUAUCGGCUAUAG
db_A: -..((((.(((......)))..))))....
Cons: .(.((((.((((....))))..)))).)..
db_B: ...((((.(((......)))..))))-...
V_B:  AAUGAUCUAUGAACUGUUAUCUGAUU-UAG
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

### Example input sequences

Example input files can be found in `src/Tests/Data`.

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
