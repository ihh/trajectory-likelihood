
# Boost
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
BOOST_PREFIX = /usr
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/program_options.hpp))
BOOST_PREFIX = /usr/local
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/program_options.hpp))
BOOST_PREFIX = /usr/local/homebrew/opt/boost
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/program_options.hpp))
BOOST_PREFIX =
endif
endif
endif

BOOST_FLAGS := -I$(BOOST_PREFIX)/include
BOOST_LIBS := -L$(BOOST_PREFIX)/lib -lboost_program_options

# Compiler
CPP = clang++

# Flags
CPP_FLAGS = -std=c++11 -g -O3 $(BOOST_FLAGS)
LD_FLAGS = $(BOOST_LIBS)

# Targets
trajeclike: $(wildcard *.cpp) $(wildcard *.h) Makefile
	$(CPP) $(CPP_FLAGS) $(LD_FLAGS) -o $@ $(wildcard *.cpp)

clean:
	rm trajeclike
