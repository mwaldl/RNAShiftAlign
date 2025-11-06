bin_PROGRAMS = rnashiftalign

AM_DEFAULT_SOURCE_EXT = .cc

rnashiftalign_SOURCES = RNAShiftAlign.cc \
	RNAShiftAlign/shift_scoring.cc

# Test programs
BINTESTS = test_shiftmatrix_m test_shift_scoring
SCRIPTTESTS =

test_shiftmatrix_m_SOURCES = Tests/test_main.cc Tests/test_shiftmatrix_m.cc
test_shift_scoring_SOURCES = Tests/test_main.cc Tests/test_shift_scoring.cc \
	RNAShiftAlign/shift_scoring.cc

TESTS = $(BINTESTS) $(SCRIPTTESTS)
check_PROGRAMS = $(BINTESTS)

EXTRA_DIST = $(SCRIPTTESTS)
