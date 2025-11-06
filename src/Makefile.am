bin_PROGRAMS = RNAShiftAlign

AM_DEFAULT_SOURCE_EXT = .cc

RNAShiftAlign_SOURCES = RNAShiftAlign.cc

# Test programs
BINTESTS = test_shiftmatrix_m
SCRIPTTESTS =

test_shiftmatrix_m_SOURCES = Tests/test_main.cc Tests/test_shiftmatrix_m.cc

TESTS = $(BINTESTS) $(SCRIPTTESTS)
check_PROGRAMS = $(BINTESTS)

EXTRA_DIST = $(SCRIPTTESTS)
