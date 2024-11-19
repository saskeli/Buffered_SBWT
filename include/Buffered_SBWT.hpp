#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/rank_support_v.hpp"
#include "kmer.hpp"
#include "throwing_streams.hpp"
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
  uint64_t buffer_limit_;

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
    bool b_pred = false;

    bool operator<(const B_elem& rhs) const { return kmer < rhs.kmer; }
  };

  static_assert(k > 0);
  static_assert(precalc_k <= k);
  static_assert(precalc_k <= 20);
  static_assert(precalc_k > 0);

  constexpr uint64_t gigs_to_belems(double gigs) const {
    gigs *= 1024 * 1024 * 1024;
    gigs /= sizeof(B_elem);
    return uint64_t(gigs) / 2;
  }

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
    if constexpr (end - start > precalc_k) {
      KM_t s_km = kmer;
      for (uint16_t i = 0; i < start; ++i) {
        s_km = s_km.suf();
      }
      auto precalc_idx = s_km.get_first_v<precalc_k>();
      auto I = kmer_prefix_precalc_[precalc_idx];
      a = I.first;
      b = I.second;
      for (uint16_t i = precalc_k; i < end - start; ++i) {
        fl(a, b, s_km.get_v(i));
        if constexpr (short_circuit) {
          if (b <= a) {
            return {0, 0};
          }
        }
      }
      return {a, b};
    }
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

  void setup_buffer() {
    new_bits_ = {sdsl::bit_vector(n_nodes_), sdsl::bit_vector(n_nodes_),
                 sdsl::bit_vector(n_nodes_), sdsl::bit_vector(n_nodes_)};
  }

  void compute_buffer_edges() {
    __gnu_parallel::sort(buffer.begin(), buffer.end());
    uint64_t trg = 0;
    for (uint64_t src = 1; src < buffer.size(); src++) {
      if (buffer[trg].kmer != buffer[src].kmer) {
        buffer[++trg] = buffer[src];
      }
    }
    ++trg;
    buffer.resize(trg);

    uint64_t a_start = 0;
    uint64_t i = 0;
    while (i < buffer.size() && buffer[i].kmer.get_v(k - 1) == 0) {
      ++i;
    }
    uint64_t c_start = i;
    while (i < buffer.size() && buffer[i].kmer.get_v(k - 1) == 1) {
      ++i;
    }
    uint64_t g_start = i;
    while (i < buffer.size() && buffer[i].kmer.get_v(k - 1) == 2) {
      ++i;
    }
    uint64_t t_start = i;
    uint64_t ai = a_start;
    uint64_t ci = c_start;
    uint64_t gi = g_start;
    uint64_t ti = t_start;
    a_start = trg;
    Kmer_t prev;
    for (i = 0; i < trg; ++i) {
      Kmer_t ip = buffer[i].kmer.suf();
      if (i == 0 || ip != prev) {
        buffer[i].group_head = true;
      } else {
        continue;
      }
      prev = ip;
      while (ai < c_start) {
        auto ec = buffer[ai].kmer.pref();
        if (ip == ec) {
          buffer[i].edge[0] = 1;
          buffer[ai].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ai;
      }
      while (ci < g_start) {
        auto ec = buffer[ci].kmer.pref();
        if (ip == ec) {
          buffer[i].edge[1] = 1;
          buffer[ci].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ci;
      }
      while (gi < t_start) {
        auto ec = buffer[gi].kmer.pref();
        if (ip == ec) {
          buffer[i].edge[2] = 1;
          buffer[gi].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++gi;
      }
      while (ti < a_start) {
        auto ec = buffer[ti].kmer.pref();
        if (ip == ec) {
          buffer[i].edge[3] = 1;
          buffer[ti].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ti;
      }
    }
  }

  void commit() {
    setup_buffer();
    compute_buffer_edges();
    Dummy_trie dummies(*this);
#pragma omp parallel for
    for (uint64_t buffer_idx = 0; buffer_idx < buffer.size(); ++buffer_idx) {
      if (buffer[buffer_idx].group_head == false) {
        continue;
      }
      // TODO: Given a suffix group bv, potential location, and predecessor
      // suffix group, can this be done better?
      auto I = search_kmer<k, 1>(buffer[buffer_idx].kmer);
      // This (I) is the suffix group in the static_sbwt.
      // If the suffix group is empty {-1, -1}, then the k-mer
      //  will be its own suffix group, unless there is another
      //  element with the same suffix in the buffer.
      // Else the element is part of an existing suffix group.
      //  depending on other buffer elements and colex ordering.
      // Note that dummies are always singleton suffix groups:
      //  If the dummy is part of a suffix group, any descendent of the
      //  dummy, could be the descendent of some other node, making the
      //  dummy superfluous, which is a contradiction, since dummies are
      //  only added if needed.
      if (I.first != I.second) {
        uint64_t a = I.first;  // Suffix group leader.
        uint64_t a_w = a / 64;
        uint64_t a_b = uint64_t(1) << (a % 64);
        uint64_t a_nb = ~a_b;
        uint64_t d_tree_idx;
#pragma omp critical
        { d_tree_idx = dummies.mark_for_removal(buffer[buffer_idx].kmer); }
        if (d_tree_idx > 0) {
          // We are replacing a dummy with a suffix group from the
          // buffer. The suffix group start from the buffer does not
          // need to change.
          for (uint16_t i = 0; i < 4; ++i) {
            buffer[buffer_idx].edge[i] |= dummies.nodes[d_tree_idx].out[i] > 0;
          }
        } else if (buffer[buffer_idx].source == a) {
          // the buffer element will become the suffix group leader.
          // Steal children from the old leader.
          uint64_t* d = suffix_group_starts_.data();
#pragma omp atomic
          d[a_w] = d[a_w] & a_nb;
          for (uint16_t ci = 0; ci < 4; ++ci) {
            if (bits_[ci][a]) {
              buffer[buffer_idx].edge[ci] = true;
              d = new_bits_[ci].data();
#pragma omp atomic
              d[a_w] = d[a_w] | a_b;
            }
          }
        } else {
          // The old suffix group leader stays.
          // Give children to old leader.
          for (uint16_t ci = 0; ci < 4; ++ci) {
            if (buffer[buffer_idx].edge[ci] > 0) {
              uint64_t* d = new_bits_[ci].data();
#pragma omp atomic
              d[a_w] = d[a_w] | a_b;
            }
          }
          buffer[buffer_idx].edge = {0, 0, 0, 0};
          buffer[buffer_idx].group_head = false;
        }
      }
      // If the buffer element has no predecessor in buffer, the incoming edge
      // needs to either be added to an existing predecessor in sbwt or a new
      // dummy path needs to be created
      if (buffer[buffer_idx].b_pred == false) {
        if (buffer[buffer_idx].pred_size > 0) {
          uint64_t a_w = a / 64;
          uint64_t a_b = ONE << (a % 64);
          uint64_t* d = new_bits_[buffer[buffer_idx].kmer.get_v(k - 1)].data();
#pragma omp atomic
          d[a_w] = d[a_w] < a_b;
        } else {
#pragma omp critical
          { dummies.add(buffer[buffer_idx].kmer, *this); }
        }
      }
    }

    std::vector<uint64_t> removables;
    std::vector<typename decltype(dummies)::Dummy_trie_node> addables;
    dummies.update_removals(removables, addables, *this);

#pragma omp parallel for
    for (uint16_t i = 0; i < 4; ++i) {
      new_bits_[0] ^= bits_[0];
    }

    sdsl::bit_vector old_starts;
    old_starts = suffix_group_starts_;

    uint64_t o_size = n_nodes_;
    uint64_t n_size =
        o_size + buffer.size() + addables.size() - removables.size();
    for (uint16_t i = 0; i < 4; ++i) {
      bits_[i].bit_resize(n_size);
    }
    suffix_group_starts_.bit_resize(n_size);

    n_nodes_ = n_size;
    n_kmers_ += buffer.size();

    uint64_t buffer_index = 0;
    uint64_t sbwt_index = 0;
    uint64_t dummy_index = 0;
    uint64_t r_dummy_index = 0;
    uint64_t write_index = 0;
    uint64_t next_i = std::min(
        {buffer.size() > 0 ? buffer[buffer_index].source : n_size,
         addables.size() > 0 ? addables[dummy_index].sbwt_index : n_size,
         removables.size() > 0 ? removables[r_dummy_index] : n_size});
    while (buffer_index < buffer.size() || sbwt_index < o_size ||
           dummy_index < addables.size()) {
      if (buffer_index < buffer.size()) {
        auto be = buffer[buffer_index];
        if (be.source == sbwt_index) {
          if (dummy_index < addables.size()) {
            auto nd = addables[dummy_index];
            if (nd.sbwt_index == sbwt_index && nd.kmer <= be.kmer) {
              for (uint16_t i = 0; i < 4; ++i) {
                bits_[i][write_index] = nd.out[i];
              }
              suffix_group_starts_[write_index] = true;
              ++dummy_index;
              ++write_index;
              next_i = std::min(
                  {buffer.size() > buffer_index ? buffer[buffer_index].source
                                                : o_size,
                   addables.size() > dummy_index
                       ? addables[dummy_index].sbwt_index
                       : o_size,
                   removables.size() > r_dummy_index ? removables[r_dummy_index]
                                                     : o_size});
              continue;
            }
          }
          for (uint16_t i = 0; i < 4; ++i) {
            bits_[i][write_index] = be.edge[i];
          }
          suffix_group_starts_[write_index] = be.group_head;
          ++buffer_index;
          ++write_index;
          next_i = std::min(
              {buffer.size() > buffer_index ? buffer[buffer_index].source
                                            : o_size,
               addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                             : o_size,
               removables.size() > r_dummy_index ? removables[r_dummy_index]
                                                 : o_size});
          continue;
        }
      }
      if (dummy_index < addables.size()) {
        auto nd = addables[dummy_index];
        if (nd.sbwt_index == sbwt_index) {
          for (uint16_t i = 0; i < 4; ++i) {
            bits_[i][write_index] = nd.out[i];
          }
          suffix_group_starts_[write_index] = true;
          ++dummy_index;
          ++write_index;
          next_i = std::min(
              {buffer.size() > buffer_index ? buffer[buffer_index].source
                                            : o_size,
               addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                             : o_size,
               removables.size() > r_dummy_index ? removables[r_dummy_index]
                                                 : o_size});
          continue;
        }
      }
      if (r_dummy_index < removables.size() &&
          removables[r_dummy_index] == sbwt_index) {
        ++r_dummy_index;
        ++sbwt_index;
        next_i = std::min(
            {buffer.size() > buffer_index ? buffer[buffer_index].source
                                          : o_size,
             addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                           : o_size,
             removables.size() > r_dummy_index ? removables[r_dummy_index]
                                               : o_size});
        continue;
      }
      uint64_t copy_count = next_i - sbwt_index;
      uint64_t w = write_index / 64;
      uint8_t w_o = write_index % 64;
      write_index += copy_count;
      uint64_t r = sbwt_index / 64;
      uint8_t r_o = sbwt_index % 64;
      sbwt_index = next_i;
      std::array<uint64_t*, 4> w_ptr = {
          bits_[0].data() + w, bits_[1].data() + w, bits_[2].data() + w,
          bits_[3].data() + w};
      uint64_t* s_w_ptr = suffix_group_starts_.data() + w;
      std::array<uint64_t*, 4> r_ptr = {
          new_bits_[0].data() + r, new_bits_[1].data() + r,
          new_bits_[2].data() + r, new_bits_[3].data() + r};
      uint64_t* s_r_ptr = old_starts.data() + r;
      while (copy_count > 64) {
        for (uint16_t i = 0; i < 4; ++i) {
          w = sdsl::bits::read_int(r_ptr[i]++, r_o);
          sdsl::bits::write_int(w_ptr[i]++, w, w_o);
        }
        w = sdsl::bits::read_int(s_r_ptr++, r_o);
        sdsl::bits::write_int(s_w_ptr++, w, w_o);
        copy_count -= 64;
      }
      for (uint16_t i = 0; i < 4; ++i) {
        w = sdsl::bits::read_int(r_ptr[i], r_o, copy_count);
        sdsl::bits::write_int(w_ptr[i], w, w_o, copy_count);
      }
      w = sdsl::bits::read_int(s_r_ptr, r_o, copy_count);
      sdsl::bits::write_int(s_w_ptr, w, w_o, copy_count);
    }

    for (uint16_t i = 0; i < 4; ++i) {
      sdsl::util::init_support(rank_supports_[i], &(bits_[i]));
    }

    for (uint16_t i = 1; i <= 4; ++i) {
      c_table_[i] = rank_supports_[i].rank(n_size) + c_table_[i - 1];
    }

    compute_k_mer_precalc();
    buffer.clear();
  }

 public:
  Buffered_SBWT(double buffer_gigs = 0.5)
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
    buffer_limit_ = gigs_to_belems(buffer_gigs);
  }

  Buffered_SBWT(const std::string& file, double buffer_gigs = 0.5)
      : Buffered_SBWT(buffer_gigs) {
    load(file);
  }

  template <class IS>
  Buffered_SBWT(IS& in, double buffer_gigs = 0.5) : Buffered_SBWT(buffer_gigs) {
    load(in);
  }

  Buffered_SBWT(const Buffered_SBWT&) = delete;
  Buffered_SBWT(Buffered_SBWT&&) = delete;
  Buffered_SBWT& operator=(Buffered_SBWT&&) = delete;
  Buffered_SBWT& operator=(Buffered_SBWT) = delete;

  uint64_t serialize(const std::string& file) const {
    throwing_ofstream out(file, ios::binary);
    return serialize(out);
  }

  template <class OS>
  uint64_t serialize(OS& out_stream) const {
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

  uint64_t serialize_old_format(const std::string& file) const {
    throwing_ofstream out(file, ios::binary);
    return serialize_old_format(out);
  }

  template<class OS>
  uint64_t serialize_old_format(OS& out_stream) const {
    uint64_t written = 0;
    written += serialize_string(SBWT_VARIANT, out_stream);
    std::string old_version = "v0.1";
    written += serialize_string(old_version, out_stream);
    for (uint16_t i = 0; i < 4; ++i) {
        bits_[i].serialize(out_stream);
    }
    for (uint16_t i = 0; i < 4; ++i) {
        rank_supports_[i].serialize(out_stream);
    }
    written += suffix_group_starts_.serialize(out_stream);
    std::vector<uint64_t> C;
    for (uint16_t i = 0; i < 4; ++i) {
        C.push_back(c_table_[i]);
    }
    written += serialize_std_vector(C, out_stream);
    std::vector<std::pair<int64_t, int64_t>> old_precalc;
    for (auto P : kmer_prefix_precalc_) {
        old_precalc.push_back(P.first == P.second ? {-1, -1}, {P.first, P.second - 1});
    }
    written += serialize_std_vector(old_precalc, out_stream);
    int64_t writable_v = precalc_k;
    out_stream.write(reinterpret_cast<char*>(&writable_v), sizeof(writable_v));
    written += sizeof(writable_v);

    writable_v = n_nodes_;
    out_stream.write(reinterpret_cast<char*>(&writable_v), sizeof(writable_v));
    written += sizeof(writable_v);

    writable_v = n_kmers_;
    out_stream.write(reinterpret_cast<char*>(&writable_v), sizeof(writable_v));
    written += sizeof(writable_v);

    writable_v = k;
    out_stream.write(reinterpret_cast<char*>(&writable_v), sizeof(writable_v));
    written += sizeof(writable_v);

    return written;
  }

  void load(const std::string& file) {
    throwing_ifstream in(file, ios::binary);
    load(in);
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
#ifdef DEBUG
    if (not is_valid()) {
      exit(1);
    }
#endif
  }

  template <class addables_t>
  uint64_t add_all(addables_t& addables) {
    uint64_t total_k_mers = 0;
    bool cont = true;
    while (cont) {
      uint64_t k_mer_count = 0;
// parallel
#pragma omp parallel
      {
        std::vector<B_elem> addable_k_mers;
        while (true) {
          bool do_break = false;
          std::string addable;
#pragma omp critical
          {
            do_break = not addables.get(addable);
            cont = cont && not do_break;
            if (not do_break) {
              total_k_mers += addable.size() - k + 1;
            }
          }
          if (do_break) {
            if (addable_k_mers.size() > 0) {
#pragma omp critical
              {
                buffer.insert(buffer.end(), addable_k_mers.begin(),
                              addable_k_mers.end());
              }
            }
            break;
          }
          add_string(addable, addable_k_mers);
#pragma omp atomic
          k_mer_count = k_mer_count + addable_k_mers.size();

          if (k_mer_count >= buffer_limit) {
            if (addable_k_mers.size() > 0) {
              std::sort(addable_k_mers.begin(), addable_k_mers.end());
              uint64_t trg = 0;
              for (uint64_t i = 0; i < addable_k_mers.size(); ++i) {
                if (addable_k_mers[trg].kmer != addable_k_mers[i].kmer) {
                  addable_k_mers[++trg] = addable_k_mers[i];
                }
              }
              ++trg;
              addable_k_mers.resize(trg);
#pragma omp critical
              {
                buffer.insert(buffer.end(), addable_k_mers.begin(),
                              addable_k_mers.end());
              }
            }
            break;
          }
        }
      }
      if (buffer.size() >= buffer_limit) {
        commit();
#ifdef DEBUG
        if (not is_valid()) {
          exit(1);
        }
#endif
      }
    }
    if (buffer.size() > 0) {
      commit();
#ifdef DEBUG
      if (not is_valid()) {
        exit(1);
      }
#endif
    }
    return total_k_mers;
  }

  uint64_t number_of_kmers() const { return n_kmers_; }

  int64_t search(const std::string& kmer) const {
    Kmer_t km(kmer);
    auto I = search_kmer(km);
    return I.first == I.second ? -1 : I.first;
  }

  std::vector<int64_t> streaming_search(const std::string& input) const {
    std::vector<int64_t> ret;
    Streaming_search ss(input, *this);
    uint64_t i;
    bool found;
    while (ss.next(i, found)) {
      ret.push_back(found ? i : -1);
    }
    return ret;
  }

  uint16_t get_precalc_k() const { return precalc_k; }

  uint16_t get_k() const { return k; }

  void print() {
    std::cout << "SBWT C-table, matrix, suffix group start and prefix "
                 "precalc intervals."
              << std::endl;
    std::cout << "a c g t    sgs" << std::endl;
    std::cout << c_table_[1] << " " << c_table_[2] << " " << c_table_[3] << " "
              << c_table_[4] << " <- C-table" << std::endl;
    for (uint64_t i = 0; i < n_nodes_; ++i) {
      std::cout << bits_[0][i] << " " << bits_[1][i] << " " << bits_[2][i]
                << " " << bits_[3][i] << "    " << suffix_group_starts_[i]
                << std::endl;
    }
    /*std::cout << "Prefix precald:\n";
    for (auto p : static_sbwt.kmer_prefix_precalc) {
        std::cout << p.first << ", " << p.second << std::endl;
    }//*/
  }

  bool is_valid() {
    uint64_t n = n_nodes_;
    bool ok = true;
    uint64_t c = bits_[0].size();
    if (c != n) {
      std::cerr << "subset count != A_bits.size, " << n << " != " << c
                << std::endl;
      ok = false;
    }
    c = bits_[1].size();
    if (c != n) {
      std::cerr << "subset count != C_bits.size, " << n << " != " << c
                << std::endl;
      ok = false;
    }
    c = bits_[2].size();
    if (c != n) {
      std::cerr << "subset count != G_bits.size, " << n << " != " << c
                << std::endl;
      ok = false;
    }
    c = bits_[3].size();
    if (c != n) {
      std::cerr << "subset count != T_bits.size, " << n << " != " << c
                << std::endl;
      ok = false;
    }
    c = rank_supports_[0].rank(n);
    c += rank_supports_[1].rank(n);
    c += rank_supports_[2].rank(n);
    c += rank_supports_[3].rank(n);
    if (c != n - 1) {
      std::cerr << "tree is not. number of edges should be " << n - 1 << ", is "
                << c << std::endl;
      ok = false;
    }
    return ok;
  }
};

template <uint16_t k, uint16_t precalc_k>
class Buffered_SBWT<k, precalc_k>::Streaming_search {
 private:
  const Buffered_SBWT& static_bsbwt_;
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
      static_bsbwt_.fl(a, b, kmer_.get_v(k - 1));
      loc_ = a;
      found_ = a < b;
    } else {
      auto precalc_idx = kmer_.get_first_v<precalc_k>();
      auto p = static_bsbwt_.kmer_prefix_precalc_[precalc_idx];
      uint64_t a = p.first;
      uint64_t b = p.second + 1;
      for (uint16_t i = precalc_k; i < k - 1; ++i) {
        static_bsbwt_.fl(a, b, kmer_.get_v(i));
      }
      pred_ = a;
      p_size_ = b - a;
      static_bsbwt_.fl(a, b, kmer_.get_v(k - 1));
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

template <uint16_t k, uint16_t precalc_k>
class Buffered_SBWT<k, precalc_k>::Dummy_trie {
 public:
  struct Dummy_trie_node {
    Kmer<k> kmer;
    std::array<uint64_t, 4> out;
    uint64_t sbwt_index;
    uint16_t kmer_length;
    bool keep;

    bool operator<(const Dummy_trie_node& rhs) const {
      if (kmer == rhs.kmer) {
        return kmer_length < rhs.kmer_length;
      }
      return kmer < rhs.kmer;
    }
  };

  std::vector<Dummy_trie_node> nodes;

 private:
  void walk(uint64_t index, uint16_t depth, Buffered_SBWT& sbwt) {
    if (depth == k) {
      for (uint16_t i = 0; i < 4; ++i) {
        nodes[index].out[i] = sbwt.bits_[i][nodes[index].sbwt_index];
      }
      return;
    }
    for (uint16_t i = 0; i < 4; ++i) {
      if (sbwt.bits_[i][nodes[index].sbwt_index]) {
        uint64_t idx = nodes.size();
        uint64_t sbwt_index =
            sbwt.rank_supports_[i].rank(nodes[index].sbwt_index) +
            sbwt.c_table_[i];
        nodes.push_back({nodes[index].kmer.next_v(i),
                         {0, 0, 0, 0},
                         sbwt_index,
                         depth,
                         true});
        nodes[index].out[i] = idx;
        walk(idx, depth + 1, sbwt);
      }
    }
  }

 public:
  Dummy_trie(Buffered_SBWT& sbwt) : nodes() {
    nodes.push_back({Kmer<k>(), {0, 0, 0, 0}, 0, 0, true});
    walk(0, 1, sbwt);
  }

  template <class Kmer_t>
  uint64_t mark_for_removal(Kmer_t& kmer) {
    uint64_t dummy_id = 0;
    for (uint16_t i = 1; i < k; ++i) {
      uint16_t v = kmer.get_v(i);
      dummy_id = nodes[dummy_id].out[v];
      if (dummy_id == 0) {
        return 0;
      }
    }
    nodes[dummy_id].keep = false;
    return dummy_id;
  }

  template <class Kmer_t>
  void add(Kmer_t& kmer, Buffered_SBWT& sbwt) {
    uint64_t id = 0;
    uint16_t i = 0;
    while (i < k) {
      uint16_t v = kmer.get_v(i);
      if (nodes[id].out[v] > 0) {
        id = nodes[id].out[v];
      } else {
        if (id == 0 || nodes[id].sbwt_index != 0) {
          sbwt.new_bits_[v][nodes[id].sbwt_index] = true;
        }
        if (i == k - 1) {
          nodes[id].out[v] = 1;
        }
        break;
      }
      ++i;
    }
    while (i < k - 1) {
      uint64_t idx = nodes.size();
      uint16_t v = kmer.get_v(i);
      Dummy_trie_node nd = {
          nodes[id].kmer.next_v(v), {0, 0, 0, 0}, 0, uint16_t(i + 1), true};
      if (i == k - 2) {
        nd.out[kmer.get_v(i + 1)] = 1;
      }
      nodes.push_back(nd);
      nodes[id].out[v] = idx;
      id = idx;
      ++i;
    }
  }

  template <class R_t, class A_t>
  bool mark(uint64_t node_index, uint16_t depth, R_t& removables, A_t& addables,
            Buffered_SBWT& sbwt, uint64_t sbwt_index) {
    if (depth < k) {
      nodes[node_index].keep = false;
      for (uint16_t i = 0; i < 4; ++i) {
        uint64_t c_id = nodes[node_index].out[i];
        if (c_id > 0) {
          uint64_t n_sbwt_index =
              sbwt.rank_supports_[i].rank(sbwt_index) + sbwt.c_table_[i];
          bool o_keep =
              mark(c_id, depth + 1, removables, addables, sbwt, n_sbwt_index);
          if (nodes[c_id].sbwt_index > 0 && not o_keep) {
            sbwt.new_bits_[i][sbwt_index] = true;
          }
          nodes[node_index].keep |= o_keep;
        }
      }
    }
    if (nodes[node_index].keep == false) {
      removables.push_back(nodes[node_index].sbwt_index);
      return false;
    } else if (nodes[node_index].sbwt_index == 0) {
      auto simplified = nodes[node_index];
      simplified.out = {simplified.out[0] > 0, simplified.out[1] > 0,
                        simplified.out[2] > 0, simplified.out[3] > 0};
      simplified.sbwt_index = sbwt_index;
      addables.push_back(simplified);
    }
    return true;
  }

  template <class R_t, class A_t>
  void update_removals(R_t& removables, A_t& addables, Buffered_SBWT& sbwt) {
    for (uint16_t i = 0; i < 4; ++i) {
      uint64_t b = nodes[0].out[i];
      if (b > 0) {
        bool keep = mark(b, 2, removables, addables, sbwt, sbwt.c_table_[i]);
        if (nodes[b].sbwt_index > 0 && not keep) {
          sbwt.new_bits_[i][0] = true;
        }
      }
    }
    std::sort(removables.begin(), removables.end());
    std::sort(addables.begin(), addables.end());
  }

  void print() {
    std::cout << "Dummy tire:" << std::endl;
    uint64_t idx = 0;
    for (auto de : nodes) {
      std::cout << (idx++) << ": " << de.kmer.to_string() << ": " << de.out[0]
                << " " << de.out[1] << " " << de.out[2] << " " << de.out[3]
                << ", " << de.sbwt_index << std::endl;
    }
  }

  template <class vec>
  bool validation_walk(uint64_t index, vec& sources, uint16_t depth) {
    if (depth == k) {
      return true;
    }
    for (uint16_t i = 0; i < 4; ++i) {
      uint64_t trg = nodes[index].out[i];
      if (trg > 0) {
        if (sources[trg] >= 0) {
          return false;
        }
        sources[trg] = index;
        if (not validation_walk(trg, sources, depth + 1)) {
          return false;
        }
      }
    }
    return true;
  }

  bool is_valid() {
    std::vector<int64_t> sources(nodes.size(), -1);
    sources[0] = 0;
    if (not validation_walk(0, sources, 1)) {
      return false;
    }
    for (auto v : sources) {
      if (v == -1) {
        return false;
      }
    }
    return true;
  }
};

}  // namespace sbwt