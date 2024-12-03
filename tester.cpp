#include <omp.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>

#include "Buffered_SBWT.hpp"

#define K 8

#define PRECALC_K 3

typedef sbwt::Buffered_SBWT<K, PRECALC_K> buf_t;

typedef std::unordered_set<std::string> ss_t;

void add_all(ss_t& ss, const std::string& s) {
  for (uint64_t i = 0; i < s.size() - K + 1; ++i) {
    ss.insert(s.substr(i, K));
  }
}

void del_all(ss_t& ss, const std::string& s) {
  for (uint64_t i = 0; i < s.size() - K + 1; ++i) {
    ss.erase(s.substr(i, K));
  }
}

int main(int argc, char const* argv[]) {
  omp_set_num_threads(1);
  uint16_t n = K + 20 - 1;
  std::mt19937 rnd(1337);
  std::string a;
  std::string b;
  std::string c;
  std::string d;
  std::array<char, 4> alpha = {'A', 'C', 'G', 'T'};

  if (argc >= 3) {
    buf_t buf;
    ss_t control;
    a = argv[1];
    b = argv[2];
    buf.add(a);
    buf.add(b);
    add_all(control, a);
    add_all(control, b);
    std::cout << "ADDED: " << std::endl;
    buf.print();
    for (auto s : control) {
      std::cout << "ss " << s << std::endl;
    }
    if (not buf.is_valid()) {
      exit(1);
    }

    for (uint16_t i = 3; i < argc; ++i) {
      c = argv[i];
      buf.del(c);
      del_all(control, c);

      std::cout << "\nDEL " << i << ": " << c << std::endl;
      buf.print();
      for (auto s : control) {
        std::cout << "ss " << s << std::endl;
      }
      if (not buf.is_valid()) {
        exit(1);
      }
    }

    std::cout << "OK?" << std::endl;
    exit(0);
  }

  for (uint16_t i = 0; i < n; ++i) {
    a.push_back(alpha[rnd() % 4]);
    b.push_back(alpha[rnd() % 4]);
  }
  uint64_t counter = 0;
  while (true) {
    ++counter;
    if (counter % 10000 == 0) {
      std::cout << (counter / 1000) << "k\r" << std::flush;
    }
    buf_t buf;
    buf.add(a);
    buf.add(b);

    ss_t control;
    add_all(control, a);
    add_all(control, b);

    ss_t d_set;
    for (std::string s : buf.reconstruct_all_kmers()) {
      if (s.find('$') == std::string::npos) {
        d_set.insert(s);
      }
    }

    if (control != d_set) {
      std::cerr << "\nFailed with: " << a << " " << b << " (" << c << ")"
                << std::endl;
      exit(1);
    }

    uint16_t len = (rnd() % (n - K + 1)) + K;
    if (rnd() % 1) {
      c = a.substr(rnd() % (n - len + 1), len);
    } else {
      c = b.substr(rnd() % (n - len + 1), len);
    }
    buf.del(c);
    del_all(control, c);
    d_set.clear();
    for (std::string s : buf.reconstruct_all_kmers()) {
      if (s.find('$') == std::string::npos) {
        d_set.insert(s);
      }
    }
    if (not buf.is_valid() || d_set != control) {
      std::cerr << "\nFailed with: " << a << " " << b << " " << c << std::endl;
      exit(1);
    }

    len = (rnd() % (n - K + 1)) + K;
    if (rnd() % 1) {
      d = a.substr(rnd() % (n - len + 1), len);
    } else {
      d = b.substr(rnd() % (n - len + 1), len);
    }
    buf.del(d);
    del_all(control, d);
    d_set.clear();
    for (std::string s : buf.reconstruct_all_kmers()) {
      if (s.find('$') == std::string::npos) {
        d_set.insert(s);
      }
    }
    if (not buf.is_valid() || d_set != control) {
      std::cerr << "\nFailed with: " << a << " " << b << " " << c << " " << d
                << std::endl;
      exit(1);
    }

    buf.del(a);
    buf.del(b);
    auto r = buf.reconstruct_all_kmers();
    if (not buf.is_valid() || r.size() != 1) {
      std::cerr << "\nEmpty fail with: " << a << " " << b << " " << c << " "
                << d << " " << a << " " << b << std::endl;
      exit(1);
    }
    if (counter % 2) {
      for (uint16_t i = 0; i < n; ++i) {
        a[i] = alpha[rnd() % 4];
      }
    } else {
      for (uint16_t i = 0; i < n; ++i) {
        b[i] = alpha[rnd() % 4];
      }
    }
  }

  return 0;
}
