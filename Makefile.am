AM_CXXFLAGS = -Wall -Werror -g -O2 -I$(top_srcdir)/src

# Yaggo variable
noinst_HEADERS = $(YAGGO_SOURCES)
BUILT_SOURCES = $(YAGGO_SOURCES)
DISTCLEANFILES = $(BUILT_SOURCES)

# Yaggo automatic rules with silencing
V_YAGGO = $(V_YAGGO_$(V))
V_YAGGO_ = $(V_YAGGO_$(AM_DEFAULT_VERBOSITY))
V_YAGGO_0 = @echo "  YAGGO " $@;
.yaggo.hpp:
	$(V_YAGGO)$(YAGGO) -o $@ $<
YAGGO_SOURCES = # Append all file to be built by yaggo

##############################
# Programs
##############################
bin_PROGRAMS = test
test_SOURCES = src/test.cc
noinst_HEADERS += src/test.hpp
YAGGO_SOURCES += src/test_cmdline.hpp

##############################
# Unit tests
##############################
TESTS = all_tests
check_PROGRAMS = all_tests
all_tests_CXXFLAGS = $(AM_CXXFLAGS) -I$(top_srcdir)/unittests
all_tests_LDFLAGS = -lpthread
all_tests_SOURCES = unittests/gtest/gtest-all.cc unittests/gtest/gtest_main.cc unittests/gtest/gtest.h

all_tests_SOURCES += unittests/test_trivial.cc
