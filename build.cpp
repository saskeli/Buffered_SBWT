#include <cstdint>
#include <iostream>
#include <string>

#include "Buffered_SBWT.hpp"
#include "IO_helper.hpp"

void help(std::string bin_l) {
  std::cout << ""
               "Do stuff with buffered sbwt structures.\
\
Usage: " << bin_l
            << " [OPTIONS] [sbwt1] [sbwt2] [file_list.txt]\
\
If two sbwt files are given without any .txt file, the two indexes are compared.\
\
If a .txt file and at least sbwt file are given, the fastas  referenced in the list file are added to index.\
\
Options:\
  sbwt1          input index 1 or output file if sbwt2 is not given.\
  sbwt2          output file, or index 2 for comparison.\
  file_list.txt  file containing paths to input fasta files, 1 per line.\
  -r             Also add reverse complement of input to index.\
  -n             Do not filter out N characters from input. (use if you know there are none.)\
  -t n           How many threads to run. Default = "
            << omp_get_max_threads() << "\
  -m X           Give a soft size limit (in giga bytes) for the buffer. Default X=0.5.\
  --old_format   Output in sbwt v0.1 instead of the default v0.2\n"
            << std::endl;
}

void comp(std::string a, std::string b) {
  buf_t bufa;
  std::string vara = sbwt::load_string(in.stream);
  bufa.load(in.stream);
  in = sbwt::throwing_ifstream(b, ios::binary);
  buf_t bufb;
  std::string varb = sbwt::load_string(in.stream);
  bufb.load(in.stream);
  if (vara != varb) {
    std::cerr << "Variant missmatch: " << vara << " != " << varb << std::endl;
    exit(1);
  }
  if (not bufa.is_valid()) {
    std::cerr << a << " validation failed" << std::endl;
    exit(1);
  }
  if (not bufb.is_valid()) {
    std::cerr << b << " validation failed" << std::endl;
    exit(1);
  }
  if (bufa.number_of_kmers() != bufb.number_of_kmers()) {
    std::cerr << "k-mer count missmatch: " << bufa.number_of_kmers()
              << " != " << bufb.number_of_kmers() << std::endl;
    exit(1);
  }
  if (bufa.number_of_subsets() != bufb.number_of_subsets()) {
    std::cerr << "node count (dummy count) missmatch: "
              << bufa.number_of_subsets() << " != " << bufb.number_of_subsets()
              << std::endl;
    exit(1);
  }
  if (bufa.getA().size() != bufb.getA().size()) {
    std::cerr << "A bit size missmatch" << std::endl;
    exit(1);
  }
  if (bufa.getA() != bufb.getA()) {
    std::cerr << "A bit content missmatch";
    exit(1);
  }
  if (bufa.getC().size() != bufb.getC().size()) {
    std::cerr << "C bit size missmatch" << std::endl;
    exit(1);
  }
  if (bufa.getC() != bufb.getC()) {
    std::cerr << "C bit content missmatch";
    exit(1);
  }
  if (bufa.getG().size() != bufb.getG().size()) {
    std::cerr << "G bit size missmatch" << std::endl;
    exit(1);
  }
  if (bufa.getG() != bufb.getG()) {
    std::cerr << "G bit content missmatch";
    exit(1);
  }
  if (bufa.getT().size() != bufb.getT().size()) {
    std::cerr << "T bit size missmatch" << std::endl;
    exit(1);
  }
  if (bufa.getT() != bufb.getT()) {
    std::cerr << "T bit content missmatch";
    exit(1);
  }
  std::cout << "All is fine?" << std::endl;
}

void add_files(std::string input_path, std::string output_path,
               std::string sbwt_path, bool rev_comp, bool filter_n,
               double buffer_gigs) {
  using std::chrono::duration_cast;
  using std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;

  auto t1 = high_resolution_clock::now();
  buf_t buf(8, buffer_gigs);
  std::string variant = "plain-matrix";
  if (sbwt_path.size() > 0) {
    sbwt::throwing_ifstream in(sbwt_path, ios::binary);
    variant = sbwt::load_string(in.stream);
    buf.load(in.stream);
  }

  auto t2 = high_resolution_clock::now();
  uint64_t o_size = buf.number_of_kmers();
  std::cout << sbwt_path << " loaded in "
            << double(duration_cast<nanoseconds>(t2 - t1).count()) / 1000000
            << " ms\n";
  std::cout << "    with " << o_size << " 31-mers" << std::endl;
  /*if (not buf.is_valid()) {
      exit(1);
  }*/

  std::vector<std::string> input_files = sbwt::readlines(input_path);

  io_container reader(input_files, rev_comp, filter_n);

  uint64_t offered_k_mers = buf.add_all(reader);

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

  assert(buf.is_valid());

  sbwt::throwing_ofstream out(output_path, ios::binary);
  sbwt::serialize_string(variant, out.stream);
  buf.serialize(out.stream);
}

int main(int argc, char const* argv[]) {
  uint16_t n = 40;
  bool rev_comp = false;
  bool filter_n = true;
  double buffer_gigs = 0.5;
  int num_threads = omp_get_max_threads();
  std::string in_sbwt = "";
  std::string in_files = "";
  std::string out_sbwt = "";

  for (size_t i = 1; i < argc; ++i) {
    if (std::strstr(argv[i], "-r")) {
      rev_comp = true;
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
        std::cerr << "At most one text file in the paramerters" << std::endl;
        help(argv[0]);
        exit(1);
      }
    } else {
      if (in_sbwt.size() == 0) {
        in_sbwt = argv[i];
      } else if (out_sbwt.size() == 0) {
        out_sbwt = argv[i];
      } else {
        std::cerr << "At most two index files in the parameters" << std::endl;
        help(argv[0]);
        exit(1);
      }
    }
  }
  if (num_threads != omp_get_max_threads()) {
    omp_set_num_threads(num_threads);
  }
  std::cout << "sbwt 1: " << in_sbwt << "\n"
            << "sbwt 2: " << out_sbwt << "\n"
            << "text file: " << in_files << "\n"
            << "extract reverse complements: " << rev_comp << "\n"
            << "filter N characters: " << filter_n << "\n"
            << "max buffer size (ish): " << buffer_gigs << " gigs\n"
            << "max threads: " << num_threads << std::endl;
  if (in_files.size() > 0) {
    if (in_sbwt.size() == 0) {
      std::cerr << "At least an output sbwt file is required" << std::endl;
      help(argv[0]);
      exit(1);
    }
    if (out_sbwt.size() == 0) {
      add_files(in_files, in_sbwt, "", rev_comp, filter_n, buffer_gigs);
    } else {
      add_files(in_files, out_sbwt, in_sbwt, rev_comp, filter_n, buffer_gigs);
    }
    exit(0);
  } else if (in_sbwt.size() > 0 && out_sbwt.size() > 0) {
    comp(in_sbwt, out_sbwt);
    exit(0);
  }
  help(argv[0]);
  return 0;
}