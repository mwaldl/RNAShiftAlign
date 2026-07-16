bin_PROGRAMS = rnashiftalign

AM_DEFAULT_SOURCE_EXT = .cc

rnashiftalign_SOURCES = RNAShiftAlign.cc \
	RNAShiftAlign/shift_scoring.cc \
	RNAShiftAlign/shift_aligner.cc

# Headers: not compiled on their own, but must be shipped in `make dist`
noinst_HEADERS = RNAShiftAlign/shift_aligner.hh \
	RNAShiftAlign/shift_scoring.hh \
	RNAShiftAlign/shiftmatrix_m.hh \
	RNAShiftAlign/shiftmatrix_d.hh \
	RNAShiftAlign/output_formats.hh


# Test programs
BINTESTS = test_shiftmatrix_m test_shift_aligner test_shiftmatrix_d
SCRIPTTESTS =

test_shiftmatrix_m_SOURCES = Tests/test_main.cc Tests/test_shiftmatrix_m.cc
# test_shift_scoring_SOURCES = Tests/test_main.cc Tests/test_shift_scoring.cc \
# 	RNAShiftAlign/shift_scoring.cc
test_shift_aligner_SOURCES = Tests/test_main.cc Tests/test_shift_aligner.cc \
	RNAShiftAlign/shift_scoring.cc \
	RNAShiftAlign/shift_aligner.cc
test_shiftmatrix_d_SOURCES = Tests/test_main.cc Tests/test_shiftmatrix_d.cc
# test_structure_alignment_SOURCES = Tests/test_main.cc Tests/test_aligner_align.cc \
#	RNAShiftAlign/shift_scoring.cc \
#	RNAShiftAlign/shift_aligner.cc

TESTS = $(BINTESTS) $(SCRIPTTESTS)
check_PROGRAMS = $(BINTESTS)

EXTRA_DIST = $(SCRIPTTESTS)
