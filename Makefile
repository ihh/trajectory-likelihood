
# GSL
GSL_PREFIX = $(shell gsl-config --prefix)

GSL_FLAGS = $(shell pkg-config --cflags gsl)
ifeq (, $(GSL_FLAGS))
GSL_FLAGS = -I$(GSL_PREFIX)/include
endif

GSL_LIBS = $(shell pkg-config --libs gsl)
ifeq (, $(GSL_LIBS))
GSL_LIBS = -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
endif

# Boost
# NB pkg-config support for Boost is lacking; see https://svn.boost.org/trac/boost/ticket/1094
BOOST_PREFIX = /usr
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/regex.h))
BOOST_PREFIX = /usr/local
ifeq (,$(wildcard $(BOOST_PREFIX)/include/boost/regex.h))
BOOST_PREFIX =
endif
endif

BOOST_PROGRAM_OPTIONS = program_options

BOOST_FLAGS := -I$(BOOST_PREFIX)/include
BOOST_LIBS := -L$(BOOST_PREFIX)/lib -lboost_regex -lboost_$(BOOST_PROGRAM_OPTIONS)

# SSL
SSL_FLAGS = -I/usr/local/opt/openssl/include
SSL_LIBS = -lssl -lcrypto -L/usr/local/opt/openssl/lib


# Compiler
CPP = clang++

CPP_FLAGS = -std=c++11 -g -O3 $(GSL_FLAGS) $(BOOST_FLAGS) $(SSL_FLAGS)
LD_FLAGS = -lz $(GSL_LIBS) $(BOOST_LIBS) $(SSL_LIBS) -lboss

# Targets
trajeclike: $(wildcard *.cpp)
	$(CPP) $(CPP_FLAGS) $(LD_FLAGS) -o $@ $^
