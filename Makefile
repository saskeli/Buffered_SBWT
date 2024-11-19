ifndef K
K = 31
endif

ifndef PRECALC_K
PRECALC_K = 8
endif

CFLAGS = -std=c++2b -Wall -Wextra -Wshadow -pedantic -march=native -DK=$(K) -DPRECALC_K=$(PRECALC_K)

PERF_FLAGS = -Ofast -DNDEBUG

DEBUG_FLAGS = -g -DDEBUG

INCLUDE = -isystem sdsl-lite/include -isystem SeqIO/include

HEADERS = include/Buffered_SBWT include/IO_helper.hpp include/kmer.hpp include/throwing_streams.hpp include/utils.hpp

.PHONY: clean all

.DEFAULT: all

%/%.hpp:

all: 

build: build.cpp $(HEADERS)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) build.cpp -o build
