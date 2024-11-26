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
  std::cout << R"(Compare 2 sbwt indexes.

Usage: )" << bin_l
            << R"( [OPTIONS] [sbwt1] [sbwt2]

Opens the two sbwt indexes in (v0.1 or v0.2 format) and compares them.

Options:
  sbwt1          First sbwt (the one whose compare function will be called).
  sbwt2          Second sbwt (the one passed to the first one).
  -h             Print help and terminate. Overrides other parameters.
  )" << std::endl;
}

typedef sbwt::Buffered_SBWT<K, PRECALC_K> buf_t;

void comp(std::string a, std::string b) {
  buf_t bufa(a);
  buf_t bufb(b);
  if (not bufa.is_valid()) {
    std::cerr << a << " validation failed" << std::endl;
  }
  if (not bufb.is_valid()) {
    std::cerr << b << " validation failed" << std::endl;
  }
  if (bufa.compare(bufb)) {
    std::cout << "All is fine?" << std::endl;
  } else {
    std::cout << "Is broken!" << std::endl;
    exit(1);
  }
}

int main(int argc, char const* argv[]) {
  std::string sbwt_a = "";
  std::string sbwt_b = "";

  for (size_t i = 1; i < size_t(argc); ++i) {
    if (std::strstr(argv[i], "-h")) {
      help(argv[0]);
      exit(0);
    }
    if (sbwt_a.size() == 0) {
      sbwt_a = argv[i];
    } else if (sbwt_b.size() == 0) {
      sbwt_b = argv[i];
    } else {
      std::cerr << "At most two sbwt files.\n"
                << " sbwt1 = " << sbwt_a << "\n"
                << " sbwt2 = " << sbwt_b << "\n"
                << " additional sbwt = " << argv[i] << std::endl;
      help(argv[0]);
      exit(1);
    }
  }

  std::cout << "sbwt1: " << (sbwt_a.size() ? sbwt_a : "N/A") << "\n"
            << "sbwt2: " << (sbwt_b.size() ? sbwt_b : "N/A") << "\n"
            << std::endl;
  if (sbwt_b.size() == 0) {
    std::cerr << "sbwt files are required." << std::endl;
    help(argv[0]);
    exit(0);
  }
  comp(sbwt_a, sbwt_b);
  return 0;
}