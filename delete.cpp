#include <cstdint>
#include <iostream>
#include <string>
#include <omp.h>

#include "Buffered_SBWT.hpp"
#include "IO_helper.hpp"

#ifndef K
#define K 31
#endif

#ifndef PRECALC_K
#define PRECALC_K 8
#endif

void help(std::string bin_l) {
  std::cout << R"(Remove elements from buffered sbwt structures.

Usage: )" << bin_l
            << R"( [OPTIONS] [sbwt1] [sbwt2] [file_list.txt]

The k-mers from the fastas referenced in the .txt file are removed from the input sbwt index. 
The created index is stored to the output sbwt index.

Input file list and output sbwt path are required arguments.
Input sbwt is not strictly required but this is a very complicated way to make an empty sbwt,
starting from an empty index.

Options:
  -i in_sbwt     Input index.
  out_sbwt       output file.
  file_list.txt  file containing paths to input fasta files, 1 per line.
  -r             Also remove reverse complement of input from index.
  -n             Do not filter out N characters from input. (use if you know there are none.)
  -t n           How many threads to run. Default = )"
            << omp_get_max_threads() << R"(
  -m X           Give a soft size limit (in giga bytes) for the buffer. Default X=0.5.
  --old_format   Output in sbwt v0.1 instead of the default v0.2.
  -h             Print help and terminate. Overrides other parameters.
  )" << std::endl;
}

typedef sbwt::Buffered_SBWT<K, PRECALC_K> buf_t;

void remove_files(std::string input_path, std::string output_path,
               std::string sbwt_path, bool rev_comp, bool filter_n,
               double buffer_gigs, bool old_output_format) {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;

  auto t1 = high_resolution_clock::now();
  buf_t buf(buffer_gigs);
  if (sbwt_path.size() > 0) {
    buf.load(sbwt_path);
  }

  auto t2 = high_resolution_clock::now();
  uint64_t o_size = buf.number_of_kmers();
  std::cout << sbwt_path << " loaded in "
            << double(duration_cast<nanoseconds>(t2 - t1).count()) / 1000000
            << " ms\n";
  std::cout << "    with " << o_size << " " << K << "-mers" << std::endl;
  /*if (not buf.is_valid()) {
      exit(1);
  }*/

  std::vector<std::string> input_files = sbwt::readlines(input_path);

  sbwt::io_container<K> reader(input_files, rev_comp, filter_n);

  uint64_t offered_k_mers = buf.del_all(reader);

  t2 = high_resolution_clock::now();
  double del_time = duration_cast<nanoseconds>(t2 - t1).count();
  del_time /= 1000000;
  uint64_t n_size = buf.number_of_kmers();
  std::cout << "Saw " << offered_k_mers << " in total" << std::endl;
  std::cout << "Added " << o_size - n_size << " k-mers in " << del_time << " ms"
            << std::endl;
  ;
  std::cout << del_time / (o_size - n_size) << "ms per added k-mer\n"
            << del_time / (offered_k_mers) << "ms per offered k-mer"
            << std::endl;

#ifdef DEBUG
  if (not buf.is_valid()) {
    std::cerr << "Validation fail!" << std::endl;
  }
#endif
  if (old_output_format) [[unlikely]] {
    buf.serialize_old_format(output_path);
  } else {
    buf.serialize(output_path);
  }
}

int main(int argc, char const* argv[]) {
  // uint16_t n = 40;
  bool rev_comp = false;
  bool filter_n = true;
  double buffer_gigs = 0.5;
  bool output_old_format = false;
  int num_threads = omp_get_max_threads();
  std::string in_sbwt = "";
  std::string in_files = "";
  std::string out_sbwt = "";

  for (size_t i = 1; i < size_t(argc); ++i) {
    if (std::strstr(argv[i], "-h")) {
      help(argv[0]);
      exit(0);
    }
    if (std::strstr(argv[i], "-r")) {
      rev_comp = true;
    } else if (std::strstr(argv[i], "-i")) {
      in_sbwt = argv[++i];
    } else if (std::strstr(argv[i], "--old_format")) {
      output_old_format = true;
    } else if (std::strstr(argv[i], "-n")) {
      filter_n = false;
    } else if (std::strstr(argv[i], "-m")) {
      buffer_gigs = std::stod(argv[++i]);
    } else if (std::strstr(argv[i], "-t")) {
      num_threads = std::stoi(argv[++i]);
    } else if (std::strstr(argv[i], ".txt")) {
      if (in_files.size() == 0) {
        in_files = argv[i];
      } else {
        std::cerr << "At most one text file in the parameters" << std::endl;
        help(argv[0]);
        exit(1);
      }
    } else {
      if (out_sbwt.size() == 0) {
        out_sbwt = argv[i];
      } else {
        std::cerr << "Exactly one output index file is required" << std::endl;
        help(argv[0]);
        exit(1);
      }
    }
  }
  if (num_threads != omp_get_max_threads()) {
    omp_set_num_threads(num_threads);
  }
  std::cout << "in sbwt: " << (in_sbwt.size() ? in_sbwt : "N/A") << "\n"
            << "out_sbwt: " << (out_sbwt.size() ? out_sbwt : "N/A") << "\n"
            << "text file: " << (in_files.size() ? in_files : "N/A") << "\n"
            << "extract reverse complements: " << rev_comp << "\n"
            << "filter N characters: " << filter_n << "\n"
            << "max buffer size (ish): " << buffer_gigs << " gigs\n"
            << "max threads: " << num_threads << "\n"
            << "output fomat: " << (output_old_format ? "v0.1\n" : "v0.2\n")
            << std::endl;
  if (in_files.size() > 0) {
    if (out_sbwt.size() == 0) {
      std::cerr << "Output sbwt file is required" << std::endl;
      help(argv[0]);
      exit(1);
    }
    remove_files(in_files, out_sbwt, in_sbwt, rev_comp, filter_n, buffer_gigs,
              output_old_format);
    exit(0);
  } else {
    std::cerr << "Input text file is required" << std::endl;
    help(argv[0]);
    exit(1);
  }
  return 0;
}