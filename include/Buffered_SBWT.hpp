#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../sdsl-lite/include/sdsl/bit_vectors.hpp"
#include "../sdsl-lite/include/sdsl/rank_support_v.hpp"
#include "kmer.hpp"
#include "utils.hpp"

namespace sbwt {

const static constexpr uint64_t ONE = 1;
const std::string SBWT_VERSION =
    "v0.2";  // Update this after breaking changes. This is serialized with the
             // index and checked when loading.
const std::string SBWT_VARIANT =
    "plain-matrix";  // Only supported variant of this version.

template <uint16_t k, uint16_t precalc_k = 8>
class Buffered_SBWT {
 private:
  std::array<sdsl::bit_vector, 4> bits_;
  std::array<sdsl::rank_support_v5<>, 4> rank_supports_;
  std::array<sdsl::bit_vector, 4> new_bits_;
  sdsl::bit_vector suffix_group_starts_;
  sdsl::bit_vector new_suffix_groups_;
  std::array<uint64_t, 5> c_table_;
  std::vector<std::pair<uint64_t, uint64_t>> kmer_prefix_precalc_;
  uint64_t n_nodes_;
  uint64_t n_kmers_;

  typedef Kmer<k> Kmer_t;

  class Dummy_trie;
  class Streaming_search;

  struct B_elem {
    Kmer_t kmer;
    uint64_t source;
    uint64_t pred;
    std::array<uint8_t, 4> edge;
    uint8_t prefix_suffix_group_size;
    bool group_head = false;

    bool operator<(const B_elem& rhs) const { return kmer < rhs.kmer; }
  };

  static_assert(k > 0);
  static_assert(precalc_k <= k);
  static_assert(precalc_k <= 20);
  static_assert(precalc_k > 0);

  void fl(uint64_t& a, uint64_t& b, uint16_t v) const {
#ifdef DEBUG
    if (a > n_nodes_ || b > n_nodes_) {
      std::cerr << "FL mapping failed\n"
                << "Size is " << n_nodes_ << ", a = " << a << ", b = " << b
                << std::endl;
      exit(1);
    }
#endif
    a = rank_supports_[v]->rank(a);
    b = rank_supports_[v]->rank(b);
    a += c_table_[v];
    b += c_table_[v];
  }

  template <class KM_t, uint16_t end = KM_t::KMER_LEN, uint16_t start = 0,
            bool short_circuit = true>
  std::pair<uint64_t, uint64_t> search_kmer(const KM_t& kmer) const {
    uint64_t a = 0;
    uint64_t b = n_nodes_;
    for (uint16_t i = start; i < end; ++i) {
      fl(a, b, kmer.get_v(i));
      if constexpr (short_circuit) {
        if (b <= a) {
          return {0, 0};
        }
      }
    }
    return {a, b};
  }

  void add_string(const std::string& addable, std::vector<B_elem>& vec) {
    Streaming_search ss(addable, *this);
    Kmer_t kmer;
    uint64_t pred_idx;
    uint16_t pred_size;
    uint64_t loc;
    bool found;
    while (ss.next(pred_idx, pred_size, loc, found, kmer)) {
      if (!found) {
        vec.push_back({kmer, loc, pred_idx, {0, 0, 0, 0}, pred_size});
      }
    }
  }

  void compute_k_mer_precalc() {
    Kmer<precalc_k> p_kmer;
    const constexpr uint64_t lim = ONE << (2 * precalc_k);
    for (uint64_t i = 0; i < lim; ++i) {
      kmer_prefix_precalc_[i] =
          search_kmer<Kmer<precalc_k>, precalc_k, 0, false>(p_kmer);
      ++p_kmer;
    }
  }

 public:
  Buffered_SBWT()
      : bits_(),
        rank_supports_(),
        new_bits_(),
        suffix_group_starts_(1),
        new_suffix_groups_(1),
        c_table_({1, 1, 1, 1, 1}),
        kmer_prefix_precalc_(ONE << (2 * precalc_k), {1, 1}),
        n_nodes_(1),
        n_kmers_(0) {
    for (uint16_t i = 0; i < 4; ++i) {
      bits_[i].resize(1);
      sdsl::util::init_support(rank_supports_[i], &bits_[i]);
    }
  }

  Buffered_SBWT(const Buffered_SBWT&) = delete;
  Buffered_SBWT(Buffered_SBWT&&) = delete;
  Buffered_SBWT& operator=(Buffered_SBWT&&) = delete;
  Buffered_SBWT& operator=(Buffered_SBWT) = delete;

  template <class OS>
  uint64_t serialize(OS& out_stream) {
    uint64_t written = 0;
    written += serialize_string(SBWT_VERSION, out_stream);
    for (uint16_t i = 0; i < 4; ++i) {
      written += bits_[i].serialize(out_stream);
      written += rank_supports_[i].serialize(out_stream);
    }
    written += suffix_group_starts_.serialize(out_stream);
    written += serialize_std_array(c_table_, out_stream);
    written += serialize_std_vector(kmer_prefix_precalc_, out_stream);
    out_stream.write(reinterpret_cast<char*>(&n_nodes_), sizeof(n_nodes_));
    written += sizeof(n_nodes_);
    out_stream.write(reinterpret_cast<char*>(&n_kmers_), sizeof(n_kmers_));
    written += sizeof(n_kmers_);
    uint16_t l_k = k;
    out_stream.write(reinterpret_cast<char*>(&l_k), sizeof(l_k));
    written += sizeof(l_k);
    return written;
  }

  template <class IS>
  void load(IS& in_stream) {
    std::string var = load_string(in_stream);
    if (var == SBWT_VERSION) {
      for (uint16_t i = 0; i < 4; ++i) {
        bits_[i].load(in_stream);
        rank_supports_[i].load(in_stream, &(bits_[i]));
      }
      suffix_group_starts_.load(in_stream);
      load_std_array(c_table_, in_stream);
      load_std_vector(kmer_prefix_precalc_, in_stream);
      is.read(reinterpret_cast<char*>(&n_nodes_), sizeof(n_nodes_));
      is.read(reinterpret_cast<char*>(&n_kmers_), sizeof(n_kmers_));
      uint16_t t_k;
      is.read(reinterpret_cast<char*>(&t_k), sizeof(t_k));
      if (t_k != k) {
        throw std::runtime_error(
            "Attempt to load sbwt with the wrong value for k")
      }
    } else if (var == SBWT_VARIANT) {
      var = load_string(in_stream);
      if (var == "v0.1") {
        for (uint16_t i = 0; i < 4; ++i) {
          bits_[i].load(in_stream);
        }
        for (uint16_t i = 0; i < 4; ++i) {
          rank_supports_[i].load(in_stream, &(bits_[i]));
        }
        suffix_group_starts_.load(in_stream);
        std::array<uint64_t, 4> t_arr;
        load_std_array(t_arr, in_stream);
        for (uint16_t i = 0; i < 4; ++i) {
          c_table_[i] = t_arr[i];
        }
        c_table_[4] = bits_[0].size();
        load_std_vector(kmer_prefix_precalc_, in_stream);
        int64_t t_p_k is.read(reinterpret_cast<char*>(&t_p_k), sizeof(t_p_k));
        is.read(reinterpret_cast<char*>(&n_nodes_), sizeof(n_nodes_));
        is.read(reinterpret_cast<char*>(&n_kmers_), sizeof(n_kmers_));
        int64_t t_k;
        is.read(reinterpret_cast<char*>(&t_k), sizeof(t_k));
        if (t_k != k) {
          throw std::runtime_error(
              "Attempt to load sbwt with the wrong value for k")
        }
        compute_k_mer_precalc();
      } else {
        throw std::runtime_error("Invalid version string or corrupted SBWT");
      }
    } else {
      throw std::runtime_error("Invalid variant/version or corrupted SBWT")
    }
  }

  void add(std::string& addable) {
    add_string(addable, buffer, false);
    if (buffer.size() == 0) {
      return;
    }
    commit();
  }
};

template <uint16_t k, uint16_t precalc_k>
class Buffered_SBWT<k, precalc_k>::Streaming_search {
 private:
  const Buffered_sbwt& static_bsbwt_;
  const std::string& searchable_;
  Kmer<k> kmer_;
  uint64_t pred_;
  uint64_t loc_;
  uint64_t index_ = k - 1;
  uint16_t p_size_ = 0;
  bool found_ = false;

 public:
  Streaming_search(const std::string& searchable, const Buffered_sbwt& bsbwt)
      : static_bsbwt_(bsbwt), searchable_(searchable), kmer_(searchable) {}

  bool next(uint64_t& idx, bool& found) {
    if (index_ >= searchable_.size()) {
      return false;
    }
    if (found_) {
      uint64_t a = loc_;
      uint64_t b = a + 1;
      for (uint16_t i = 0; i < 4; ++i) {
        if (static_bsbwt_.suffix_group_starts_[a]) {
          break;
        }
        --a;
      }
      uint64_t steps = static_bsbwt_.suffix_group_starts_.size() - b;
      steps = steps > 3 ? 3 : steps;
      for (uint16_t i = 0; i < steps; ++i) {
        if (static_bsbwt_.suffix_group_starts_[b]) {
          break;
        }
        ++b;
      }
      pred_ = a;
      p_size_ = b - a;
      // std::cout << "Searching " << kmer_.to_string() << std::endl;
      // std::cout << "Initial interval: " << a << ", " << b << std::endl;
      // std::cout << "Transition with: " << kmer_.get_v(k - 1) << std::endl;
      static_bsbwt_.fl(a, b, kmer_.get_v(k - 1));
      // std::cout << "Final interval: " << a << ", " << b << std::endl;
      loc_ = a;
      found_ = a < b;
    } else {
      auto precalc_idx = kmer_.get_first_v(precalc_k);
      // std::cout << "Searching " << kmer_.to_string() << std::endl;
      // std::cout << "precalc idx = " << precalc_idx << std::endl;
      auto p = static_bsbwt_.kmer_prefix_precalc_[precalc_idx];
      uint64_t a = p.first;
      uint64_t b = p.second + 1;
      // std::cout << a << ", " << b << std::endl;
      for (uint16_t i = precalc_k; i < k - 1; ++i) {
        static_bsbwt_.fl(a, b, kmer_.get_v(i));
        // std::cout << " " << a << ", " << b << " with " << kmer_.get_v(i)
        //                   << std::endl;
      }
      pred_ = a;
      p_size_ = b - a;
      static_bsbwt_.fl(a, b, kmer_.get_v(k - 1));
      // std::cout << " " << a << ", " << b << " with¤ " << kmer_.get_v(k - 1)
      // << std::endl;
      loc_ = a;
      found_ = a < b;
    }
    idx = loc_;
    found = found_;
    ++index_;
    if (index_ < searchable_.size()) [[likely]] {
      kmer_ = kmer_.next(searchable_[index_]);
    }
    return true;
  }

  bool next(uint64_t& pred, uint16_t& pred_size, uint64_t& idx, bool& found,
            Kkmer<k>& km) {
    km = kmer_;
    bool r = next(idx, found);
    pred = pred_;
    pred_size = p_size_;
    return r;
  }
};

}  // namespace sbwt