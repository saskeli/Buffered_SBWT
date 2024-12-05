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
  
Search sbwt for k-mers from fasta files.
  
Options:
  -i sbwt        Index to search.
  files          Fasta file or file containing paths to input fasta files, 1 per line.
                 Fasta can be gzipped or not. If file type is txt, read as list.
  -n             Do not filter out N characters from input. 
                 (use if you know there are none.)
  -d             Directory containing input files. The given path will be prepended to the
                 input fasta files. Use in case the input list contains only file names.
  -h             Print help and terminate.

Example: )" << bin_l
            << R"( -d data/ fof.txt -i out/fof.sbwt
)" << std::endl;
}

void search(std::string input_path, std::string sbwt_path, bool filter_n, std::string data_dir) {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;

  auto t1 = high_resolution_clock::now();
  buf_t buf;
  if (sbwt_path.size() > 0) {
    buf.load(sbwt_path);
  }

  auto t2 = high_resolution_clock::now();
  uint64_t o_size = buf.number_of_kmers();
  std::cout << sbwt_path << " loaded in "
            << double(duration_cast<nanoseconds>(t2 - t1).count()) / 1000000
            << " ms\n";
  std::cout << "    with " << o_size << " " << K << "-mers" << std::endl;

  std::vector<std::string> input_files;
  if (input_path.ends_with(".txt")) {
    input_files = sbwt::readlines(input_path, data_dir);
  } else {
    input_files.push_back(data_dir + input_path);
  }
  sbwt::ensure_exists(input_files);
  t2 = high_resolution_clock::now();

  std::string q_s;
  uint64_t hits = 0;
  uint64_t count = 0;
  if (input_files[0].ends_with(".gz")) {
    sbwt::io_container<K, true> reader(input_files, false, filter_n);

    while (reader.get(q_s)) {
      buf.streaming_search(q_s, count, hits);
    }
  } else {
    sbwt::io_container<K> reader(input_files, false, filter_n);

    while (reader.get(q_s)) {
      buf.streaming_search(q_s, count, hits);
    }
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
  std::string data_dir = "";

  for (size_t i = 1; i < size_t(argc); ++i) {
    std::string arg(argv[i]);
    if (arg == "-h") {
      help(argv[0]);
      exit(0);
    }
    if (arg == "-n") {
      filter_n = false;
    }
    if (arg == "-i") {
      in_sbwt = argv[++i];
    } else if (arg == "-d") {
      data_dir = argv[++i];
    } else {
      if (in_files.size() == 0) {
        in_files = argv[i];
      } else {
        std::cerr << "At most one text file in the paramerters\n"
                  << "current: " << in_files << "\n"
                  << "new: " << argv[i] << std::endl;
        help(argv[0]);
        exit(1);
      }
    }
  }
  std::cout << "sbwt: " << (in_sbwt.size() > 0 ? in_sbwt : "N/A") << "\n"
            << "text file: " << in_files << "\n"
            << "fasta dir: " << (data_dir.size() ? data_dir : "N/A") << "\n"
            << "filter N characters: " << filter_n << std::endl;
  if (in_files.size() > 0) {
    search(in_files, in_sbwt, filter_n, data_dir);
    exit(0);
  } else {
    std::cerr << "File(s) to search is required" << std::endl;
  }
  help(argv[0]);
  return 0;
}
