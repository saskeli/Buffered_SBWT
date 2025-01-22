#include <omp.h>

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

void help(std::string bin_l) {
  std::cout << R"(Build buffered sbwt structures.

Usage: )" << bin_l
            << R"( [OPTIONS] [sbwt1] [sbwt2] [file_list.txt]

The k-mers from the fastas referenced in the .txt file are added to the optional input sbwt index. 
The created index is stored to the output sbwt index.

Input file list and output sbwt path are required arguments.

Options:
  -i in_sbwt     Input index.
  out_sbwt       Output file.
  -f input_fasta Input fasta file or file containing paths to input fasta files, 1 per line.
                 Fasta can be gzipped or not. If file type is txt, read as list.
  -r             Also add reverse complement of input to index.
  -n             Do not filter out N characters from input. (use if you know there are none.)
  -t n           How many threads to run. Default = )"
            << omp_get_max_threads() << R"(
  -m X           Give a soft size limit (in giga bytes) for the buffer. Default X=0.5.
  --old_format   Output in sbwt v0.1 instead of the default v0.2.
  -d             Directory containing input files. The given path will be prepended to the
                 input fasta files. Use in case the input list contains only file names.
  -h             Print help and terminate. Overrides other parameters.

Example: )" << bin_l
            << R"( -r -m 8 -d data/ -f fof.txt out/fof.sbwt
)" << std::endl;
}

typedef sbwt::Buffered_SBWT<K, PRECALC_K> buf_t;

void add_files(std::string input_path, std::string output_path,
               std::string sbwt_path, bool rev_comp, bool filter_n,
               double buffer_gigs, bool old_output_format,
               std::string data_folder) {
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
  std::vector<std::string> input_files;
  if (input_path.ends_with(".txt")) {
    input_files = sbwt::readlines(input_path, data_folder);
  } else {
    input_files.push_back(data_folder + input_path);
  }
  sbwt::ensure_exists(input_files);
  uint64_t offered_k_mers;
  if (input_files[0].ends_with(".gz")) {
    sbwt::io_container<K, true> reader(input_files, rev_comp, filter_n);

    offered_k_mers = buf.add_all(reader);
  } else {
    sbwt::io_container<K> reader(input_files, rev_comp, filter_n);

    offered_k_mers = buf.add_all(reader);
  }

  t2 = high_resolution_clock::now();
  double add_time = duration_cast<nanoseconds>(t2 - t1).count();
  add_time /= 1000000;
  uint64_t n_size = buf.number_of_kmers();
  std::cout << "Saw " << offered_k_mers << " in total" << std::endl;
  std::cout << "Added " << n_size - o_size << " k-mers in " << add_time << " ms"
            << std::endl;
  ;
  std::cout << add_time / (n_size - o_size) << "ms per added k-mer\n"
            << add_time / (offered_k_mers) << "ms per offered k-mer"
            << std::endl;

#ifdef DEBUG
  if (not buf.is_valid()) {
    std::cerr << "Validation fail!" << std::endl;
  }
#endif
  t1 = high_resolution_clock::now();
  if (old_output_format) [[unlikely]] {
    buf.serialize_old_format(output_path);
  } else {
    buf.serialize(output_path);
  }
  t2 = high_resolution_clock::now();
  double ser_time = duration_cast<nanoseconds>(t2 - t1).count();
  ser_time /= 1000000;
  std::cout << "Wrote " << output_path << " with " << n_size << " " << K
            << "-mers to file in " << ser_time << " ms" << std::endl;
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
  std::string data_folder = "";

  for (size_t i = 1; i < size_t(argc); ++i) {
    std::string arg(argv[i]);
    if (arg == "-h") {
      help(argv[0]);
      exit(0);
    }
    if (arg == "-r") {
      rev_comp = true;
    } else if (arg == "-i") {
      in_sbwt = argv[++i];
    } else if (arg == "--old_format") {
      output_old_format = true;
    } else if (arg == "-n") {
      filter_n = false;
    } else if (arg == "-m") {
      buffer_gigs = std::stod(argv[++i]);
    } else if (arg == "-t") {
      num_threads = std::stoi(argv[++i]);
    } else if (arg == "-f") {
      in_files = argv[++i];
    } else if (arg == "-d") {
      data_folder = argv[++i];
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
  if (data_folder.size() > 0 && not data_folder.ends_with("/")) {
    data_folder.append("/");
  }
  std::cout << "in sbwt: " << (in_sbwt.size() ? in_sbwt : "N/A") << "\n"
            << "out_sbwt: " << (out_sbwt.size() ? out_sbwt : "N/A") << "\n"
            << "text file: " << (in_files.size() ? in_files : "N/A") << "\n"
            << "fasta dir: " << (data_folder.size() ? data_folder : "N/A")
            << "\n"
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
    add_files(in_files, out_sbwt, in_sbwt, rev_comp, filter_n, buffer_gigs,
              output_old_format, data_folder);
    exit(0);
  } else {
    std::cerr << "Input text file is required" << std::endl;
    help(argv[0]);
    exit(1);
  }
  return 0;
}