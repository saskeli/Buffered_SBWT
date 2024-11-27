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

.PHONY: clean

.DEFAULT: all

%/%.hpp:

all: build compare search delete

debug: debug_build debug_compare debug_search debug_delete

build: build.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) build.cpp -o build $(LIBS)

compare: compare.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) compare.cpp -o compare $(LIBS)

search: search.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) search.cpp -o search $(LIBS)

delete: delete.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) delete.cpp -o delete $(LIBS)

debug_build: build.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) build.cpp -o debug_build $(LIBS)

debug_compare: compare.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) compare.cpp -o debug_compare $(LIBS)

debug_search: search.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) search.cpp -o debug_search $(LIBS)

debug_delete: delete.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) delete.cpp -o debug_delete $(LIBS)

$(SDSL_A):
	$(cd sdsl-lite && cmake CMakelists.txt && make)

clean:
	rm -f build compare search delete debug_build debug_compare debug_search debug_delete