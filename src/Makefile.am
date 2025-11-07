bin_PROGRAMS = rnashiftalign

AM_DEFAULT_SOURCE_EXT = .cc

rnashiftalign_SOURCES = RNAShiftAlign.cc \
	RNAShiftAlign/shift_scoring.cc \
	RNAShiftAlign/shift_aligner.cc

# Test programs
BINTESTS = test_shiftmatrix_m test_shift_scoring test_shift_aligner
SCRIPTTESTS =

test_shiftmatrix_m_SOURCES = Tests/test_main.cc Tests/test_shiftmatrix_m.cc
test_shift_scoring_SOURCES = Tests/test_main.cc Tests/test_shift_scoring.cc \
	RNAShiftAlign/shift_scoring.cc
test_shift_aligner_SOURCES = Tests/test_main.cc Tests/test_shift_aligner.cc \
	RNAShiftAlign/shift_scoring.cc \
	RNAShiftAlign/shift_aligner.cc

TESTS = $(BINTESTS) $(SCRIPTTESTS)
check_PROGRAMS = $(BINTESTS)

EXTRA_DIST = $(SCRIPTTESTS)
