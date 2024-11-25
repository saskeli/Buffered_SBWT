ifndef K
K = 31
endif

ifndef PRECALC_K
PRECALC_K = 8
endif

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DK=$(K) -DPRECALC_K=$(PRECALC_K)  -fopenmp

PERF_FLAGS = -Ofast -DNDEBUG

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE = -I include -isystem sdsl-lite/include -isystem SeqIO/include

LIBS = -L sdsl-lite/lib -lsdsl

HEADERS = include/Buffered_SBWT.hpp include/IO_helper.hpp include/kmer.hpp include/throwing_streams.hpp include/utils.hpp

SDSL_A = sdsl-lite/lib/libsdsl.a

.PHONY: clean all

.DEFAULT: all

%/%.hpp:

all: build

build: build.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) build.cpp -o build $(LIBS)

debug_build: build.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) build.cpp -o debug_build $(LIBS)

$(SDSL_A):
	$(cd sdsl-lite && cmake CMakelists.txt && make)