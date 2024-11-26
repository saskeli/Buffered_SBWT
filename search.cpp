#include <cstdint>
#include <iostream>
#include <string>

#include "Buffered_SBWT.hpp"
#include "IO_helper.hpp"

#ifndef K
#define K 31
#endif

#ifndef PRECALC_K
#define PRECALC_K 8
#endif

typedef sbwt::Buffered_SBWT<K, PRECALC_K> buf_t;

void help(std::string bin_l) {
  std::cout << R"(Search from sbwt.
  
Usage: )" << bin_l
            << R"( [OPTIONS] [sbwt] [file_list.txt]
  
Search sbwt for k-mers from fasta files referenced in file_list.txt.
  
Options:
  sbwt           index to search.
  file_list.txt  file containing paths to input fasta files, 1 per line.
  -n             Do not filter out N characters from input. 
                 (use if you know there are none.)
  -h             Print help and terminate.)"
            << std::endl;
}

void search(std::string input_path, std::string sbwt_path, bool filter_n) {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;

  auto t1 = high_resolution_clock::now();
  buf_t buf(sbwt_path, 1);

  auto t2 = high_resolution_clock::now();
  uint64_t o_size = buf.number_of_kmers();
  std::cout << sbwt_path << " loaded in "
            << double(duration_cast<nanoseconds>(t2 - t1).count()) / 1000000
            << " ms\n";
  std::cout << "    with " << o_size << " " << K << "-mers" << std::endl;

  std::vector<std::string> input_files = sbwt::readlines(input_path);

  sbwt::io_container<K> reader(input_files, false, filter_n);

  std::string q_s;
  uint64_t hits = 0;
  uint64_t count = 0;
  while (reader.get(q_s)) {
    buf.streaming_search(q_s, count, hits);
  }

  t1 = high_resolution_clock::now();
  double q_time = duration_cast<nanoseconds>(t1 - t2).count();
  std::cout << "Saw " << count << " k-mers in total\n"
            << "matched " << hits << " k-mers in " << q_time / 1000000
            << " ms\n"
            << "Hit rate = " << double(hits) / count << "\n"
            << q_time / count << " ns per k-mer" << std::endl;
}

int main(int argc, char const* argv[]) {
  bool filter_n = true;
  std::string in_sbwt = "";
  std::string in_files = "";

  for (size_t i = 1; i < size_t(argc); ++i) {
    if (std::strstr(argv[i], "-h")) {
      help(argv[0]);
      exit(0);
    }
    if (std::strstr(argv[i], "-n")) {
      filter_n = false;
    } else if (std::strstr(argv[i], ".txt")) {
      if (in_files.size() == 0) {
        in_files = argv[i];
      } else {
        std::cerr << "At most one text file in the paramerters" << std::endl;
        help(argv[0]);
        exit(1);
      }
    } else {
      if (in_sbwt.size() == 0) {
        in_sbwt = argv[i];
      } else {
        std::cerr << "At most one index file in the parameters" << std::endl;
        help(argv[0]);
        exit(1);
      }
    }
  }
  std::cout << "sbwt: " << in_sbwt << "\n"
            << "text file: " << in_files << "\n"
            << "filter N characters: " << filter_n << std::endl;
  if (in_files.size() > 0) {
    if (in_sbwt.size() == 0) {
      std::cerr << "sbwt file is required" << std::endl;
      help(argv[0]);
      exit(1);
    }
    search(in_files, in_sbwt, filter_n);
    exit(0);
  } else {
    std::cerr << "txt file is required" << std::endl;
  }
  help(argv[0]);
  return 0;
}
