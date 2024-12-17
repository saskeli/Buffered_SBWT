#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <parallel/algorithm>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "kmer.hpp"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/rank_support_v.hpp"
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
  typedef Kmer<k> Kmer_t;

  class Dummy_trie;
  template <bool short_circuit>
  class Streaming_search;

  struct B_elem {
    Kmer_t kmer;
    uint64_t source;
    std::array<uint8_t, 4> edge;
    bool group_head = false;
    bool b_pred = false;

    bool operator<(const B_elem& rhs) const { return kmer < rhs.kmer; }

    std::string to_string() {
      std::string sb = kmer.to_string();
      sb.append(" ");
      for (auto e : edge) {
        sb.append(e ? "1" : "0").append(" ");
      }
      sb.append(std::to_string(source));
      sb.append("    ").append(std::to_string(group_head));
      sb.append(", ").append(std::to_string(b_pred));
      return sb;
    }
  };

  static_assert(k > 0);
  static_assert(precalc_k <= k);
  static_assert(precalc_k <= 20);
  static_assert(precalc_k > 0);

  std::array<sdsl::bit_vector, 4> bits_;
  std::array<sdsl::rank_support_v5<>, 4> rank_supports_;
  std::array<sdsl::bit_vector, 4> new_bits_;
  sdsl::bit_vector suffix_group_starts_;
  std::array<uint64_t, 4> c_table_;
  std::vector<std::pair<uint64_t, uint64_t>> kmer_prefix_precalc_;
  std::vector<B_elem> buffer_;
  uint64_t n_nodes_;
  uint64_t n_kmers_;
  uint64_t buffer_limit_;

  template <class elem_t>
  constexpr uint64_t gigs_to_elems(double gigs) const {
    gigs *= 1024 * 1024 * 1024;
    gigs /= sizeof(elem_t);
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
    a = rank_supports_[v].rank(a);
    b = rank_supports_[v].rank(b);
    a += c_table_[v];
    b += c_table_[v];
  }

  template <class KM_t, uint16_t end = KM_t::KMER_LEN, uint16_t start = 0,
            bool short_circuit = true>
  std::pair<uint64_t, uint64_t> search_kmer(const KM_t& kmer) const {
    uint64_t a = 0;
    uint64_t b = n_nodes_;
    if constexpr (end - start > precalc_k) {
      auto precalc_idx = kmer.template get_first_v<precalc_k + start>();
      precalc_idx >>= 2 * start;
      auto I = kmer_prefix_precalc_[precalc_idx];
      a = I.first;
      b = I.second;
      for (uint16_t i = precalc_k + start; i < end; ++i) {
        fl(a, b, kmer.get_v(i));
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
    Streaming_search<false> ss(addable, *this);
    Kmer_t kmer;
    uint64_t loc;
    bool found;
    while (ss.next(loc, found, kmer)) {
      if (!found) {
        vec.push_back({kmer, loc, {0, 0, 0, 0}});
      }
    }
  }

  void delete_string(const std::string& removable, std::vector<B_elem>& vec) {
    Streaming_search<true> ss(removable, *this);
    Kmer_t kmer;
    uint64_t loc;
    bool found;
    while (ss.next(loc, found, kmer)) {
      if (found) {
        vec.push_back({kmer, loc, {0, 0, 0, 0}});
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
    __gnu_parallel::sort(buffer_.begin(), buffer_.end());
    uint64_t trg = 0;
    for (uint64_t src = 1; src < buffer_.size(); src++) {
      if (buffer_[trg].kmer != buffer_[src].kmer) {
        buffer_[++trg] = buffer_[src];
      }
    }
    ++trg;
    buffer_.resize(trg);
#ifdef DEBUG
    std::cout << trg << " elements left after sort and deduplication."
              << std::endl;
#endif

    uint64_t a_start = 0;
    uint64_t i = 0;
    while (i < buffer_.size() && buffer_[i].kmer.get_v(k - 1) == 0) {
      ++i;
    }
    uint64_t c_start = i;
    while (i < buffer_.size() && buffer_[i].kmer.get_v(k - 1) == 1) {
      ++i;
    }
    uint64_t g_start = i;
    while (i < buffer_.size() && buffer_[i].kmer.get_v(k - 1) == 2) {
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
      Kmer_t ip = buffer_[i].kmer.suf();
      if (i == 0 || ip != prev) {
        buffer_[i].group_head = true;
      } else {
        continue;
      }
      prev = ip;
      while (ai < c_start) {
        auto ec = buffer_[ai].kmer.pref();
        if (ip == ec) {
          buffer_[i].edge[0] = 1;
          buffer_[ai].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ai;
      }
      while (ci < g_start) {
        auto ec = buffer_[ci].kmer.pref();
        if (ip == ec) {
          buffer_[i].edge[1] = 1;
          buffer_[ci].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ci;
      }
      while (gi < t_start) {
        auto ec = buffer_[gi].kmer.pref();
        if (ip == ec) {
          buffer_[i].edge[2] = 1;
          buffer_[gi].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++gi;
      }
      while (ti < a_start) {
        auto ec = buffer_[ti].kmer.pref();
        if (ip == ec) {
          buffer_[i].edge[3] = 1;
          buffer_[ti].b_pred = true;
        } else if (ec > ip) {
          break;
        }
        ++ti;
      }
    }
  }

  void commit() {
#ifdef DEBUG
    std::cerr << "Commit " << n_nodes_ << " with " << buffer_.size()
              << " new k-mers" << std::endl;
#endif
    setup_buffer();
    compute_buffer_edges();
    Dummy_trie dummies(*this);
#ifdef VERBOSE
    uint64_t v_last = 0;
#endif
#pragma omp parallel for
    for (uint64_t buffer_idx = 0; buffer_idx < buffer_.size(); ++buffer_idx) {
      // TODO: Given a suffix group bv, potential location, and predecessor
      // suffix group, can this be done better?
      auto I = search_kmer<Kmer_t, k, 1, false>(buffer_[buffer_idx].kmer);
      // This (I) is the suffix group in the static_sbwt.
      // If the suffix group is empty (I.first == I.second), then the k-mer
      //  will be its own suffix group, unless there is another
      //  element with the same suffix in the buffer.
      // Else the element is part of an existing suffix group.
      //  depending on other buffer elements and colex ordering.
      // Note that dummies are always singleton suffix groups:
      //  If the dummy is part of a suffix group, any descendent of the
      //  dummy, could be the descendent of some other node, making the
      //  dummy superfluous, which is a contradiction, since dummies are
      //  only added if needed.
      if (I.first != I.second && buffer_[buffer_idx].group_head) {
        uint64_t a = I.first;  // Suffix group leader.
        uint64_t a_w = a / 64;
        uint64_t a_b = uint64_t(1) << (a % 64);
        uint64_t a_nb = ~a_b;
        uint64_t d_tree_idx;
#pragma omp critical
        { d_tree_idx = dummies.mark_for_removal(buffer_[buffer_idx].kmer); }
        if (d_tree_idx > 0) {
          // We are replacing a dummy with a suffix group from the
          // buffer. The suffix group start from the buffer does not
          // need to change.
          for (uint16_t i = 0; i < 4; ++i) {
            buffer_[buffer_idx].edge[i] |= dummies.nodes[d_tree_idx].out[i] > 0;
          }
        } else if (buffer_[buffer_idx].source == a) {
          // the buffer element will become the suffix group leader.
          // Steal children from the old leader.
          uint64_t* d = suffix_group_starts_.data();
#pragma omp atomic
          d[a_w] = d[a_w] & a_nb;
          for (uint16_t ci = 0; ci < 4; ++ci) {
            if (bits_[ci][a]) {
              buffer_[buffer_idx].edge[ci] = true;
              d = new_bits_[ci].data();
#pragma omp atomic
              d[a_w] = d[a_w] | a_b;
            }
          }
        } else {
          // The old suffix group leader stays.
          // Give children to old leader.
          for (uint16_t ci = 0; ci < 4; ++ci) {
            if (buffer_[buffer_idx].edge[ci] > 0) {
              uint64_t* d = new_bits_[ci].data();
#pragma omp atomic
              d[a_w] = d[a_w] | a_b;
            }
          }
          buffer_[buffer_idx].edge = {0, 0, 0, 0};
          buffer_[buffer_idx].group_head = false;
        }
      }
      // If the buffer element has no predecessor in buffer, the incoming edge
      // needs to either be added to an existing predecessor in sbwt or a new
      // dummy path needs to be created
      if (buffer_[buffer_idx].b_pred == false) {
        uint64_t a = 0;
        uint64_t b = n_nodes_;
        for (uint16_t i = 0; i < k - 1; ++i) {
          fl(a, b, buffer_[buffer_idx].kmer.get_v(i));
          if (a >= b) {
            break;
          }
        }
        uint16_t pred_size = b - a;
        if (pred_size > 0) {
          uint64_t w_index = a;
          uint64_t a_w = w_index / 64;
          uint64_t a_b = ONE << (w_index % 64);
          uint16_t transition_char = buffer_[buffer_idx].kmer.get_v(k - 1);
          uint64_t* d = new_bits_[transition_char].data();
#pragma omp atomic
          d[a_w] = d[a_w] | a_b;
        } else {
#pragma omp critical
          { dummies.add(buffer_[buffer_idx].kmer, *this); }
        }
      }
#ifdef VERBOSE
      for (uint64_t p_i = v_last; p_i <= buffer_idx; ++p_i) {
        std::cerr << p_i << ": " << buffer_[p_i].to_string() << std::endl;
      }
      v_last = buffer_idx + 1;
#endif
    }

    std::vector<uint64_t> removables;
    std::vector<typename decltype(dummies)::Dummy_trie_node> addables;
    dummies.update_removals(removables, addables, *this);

#ifdef VERBOSE
    std::cerr << "Removable dummy indexes: " << std::endl;
    for (auto rd : removables) {
      std::cerr << rd << std::endl;
    }

    std::cerr << "Addable dummies: " << std::endl;
    for (auto ad : addables) {
      std::cerr << ad.to_string() << std::endl;
    }
#endif

#pragma omp parallel for
    for (uint16_t i = 0; i < 4; ++i) {
      new_bits_[i] ^= bits_[i];
    }

    sdsl::bit_vector old_starts = suffix_group_starts_;

    uint64_t o_size = n_nodes_;
    uint64_t n_size =
        o_size + buffer_.size() + addables.size() - removables.size();
    for (uint16_t i = 0; i < 4; ++i) {
      bits_[i].bit_resize(n_size);
    }
    suffix_group_starts_.bit_resize(n_size);

#ifdef DEBUG
    uint64_t total_compute = 0;
    for (uint16_t ci = 0; ci < 4; ++ci) {
      for (uint64_t i = 0; i < n_nodes_; ++i) {
        total_compute += new_bits_[ci][i];
      }
      for (auto be : buffer_) {
        total_compute += be.edge[ci] > 0;
      }
      for (auto de : addables) {
        total_compute += de.out[ci] > 0;
      }
      for (auto de : removables) {
        total_compute -= new_bits_[ci][de];
      }
    }
    if (total_compute != n_size - 1) {
      std::cerr << "Pre-merge edge count missmatch, should be " << n_size - 1
                << " is " << total_compute << std::endl;
    }
#endif

    n_nodes_ = n_size;
    n_kmers_ += buffer_.size();

    uint64_t buffer_index = 0;
    uint64_t sbwt_index = 0;
    uint64_t dummy_index = 0;
    uint64_t r_dummy_index = 0;
    uint64_t write_index = 0;
    uint64_t next_i = std::min(
        {buffer_.size() > 0 ? buffer_[buffer_index].source : n_size,
         addables.size() > 0 ? addables[dummy_index].sbwt_index : n_size,
         removables.size() > 0 ? removables[r_dummy_index] : n_size});
    while (buffer_index < buffer_.size() || sbwt_index < o_size ||
           dummy_index < addables.size()) {
      if (buffer_index < buffer_.size()) {
        auto be = buffer_[buffer_index];
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
                  {buffer_.size() > buffer_index ? buffer_[buffer_index].source
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
              {buffer_.size() > buffer_index ? buffer_[buffer_index].source
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
              {buffer_.size() > buffer_index ? buffer_[buffer_index].source
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
            {buffer_.size() > buffer_index ? buffer_[buffer_index].source
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

    for (uint16_t i = 1; i < 4; ++i) {
      c_table_[i] = rank_supports_[i - 1].rank(n_size);
      c_table_[i] += c_table_[i - 1];
    }

    compute_k_mer_precalc();
    buffer_.clear();
    for (uint16_t i = 0; i < 4; ++i) {
      new_bits_[i].bit_resize(0);
    }
#ifdef DEBUG
    if (not is_valid()) {
      std::cerr << "Validation fail!" << std::endl;
    }
#endif
  }

  void r_commit() {
#ifdef DEBUG
    std::cerr << "Delete at " << n_nodes_ << " with " << buffer_.size()
              << " removable k-mers" << std::endl;
#endif
#ifdef VERBOSE
    print();
#endif
    setup_buffer();
    compute_buffer_edges();
    Dummy_trie dummies(*this);
#pragma omp parallel for
    for (uint64_t buffer_idx = 0; buffer_idx < buffer_.size(); ++buffer_idx) {
// Make sure no surviving child gets orphaned.
#ifdef VERBOSE
      std::cerr << buffer_[buffer_idx].to_string() << std::endl;
#endif
      if (buffer_[buffer_idx].group_head) {
        if (not suffix_group_starts_[buffer_[buffer_idx].source]) [[unlikely]] {
          uint64_t sgl = buffer_[buffer_idx].source - 1;
          for (uint64_t i = 0; i < 3; ++i) {
            if (suffix_group_starts_[sgl]) [[likely]] {
              break;
            }
            sgl--;
          }
          uint64_t w_i = sgl / 64;
          uint64_t w_b = ONE << (sgl % 64);
          for (uint16_t i = 0; i < 4; ++i) {
            if (buffer_[buffer_idx].edge[i]) {
              uint64_t* d = new_bits_[i].data();
#pragma omp atomic
              d[w_i] = d[w_i] | w_b;
            }
          }
        } else {
          std::array<bool, 4> edges = {false, false, false, false};
          bool replace = false;
          for (uint16_t i = 0; i < 4; ++i) {
            edges[i] = bits_[i][buffer_[buffer_idx].source] !=
                       buffer_[buffer_idx].edge[i];
            replace |= edges[i];
          }
          // if replace is false everything will be fine with naive removal.
          if (replace) [[unlikely]] {
            // We need to put at least one child somewhere else in this suffix
            // group.
            uint64_t trg = buffer_[buffer_idx].source + 1;
            uint64_t sib = buffer_idx + 1;
            for (uint64_t i = 1; i < 4; ++i) {
              if (trg >= n_nodes_) [[unlikely]] {
                // There is nowhere more to go, since we are out of sbwt.
                break;
              }
              if (suffix_group_starts_[trg]) [[likely]] {
                // There is nowhere more to go, since we are out of suffix
                // group.
                break;
              }
              // There is another element in the suffix group.
              if (sib < buffer_.size() && buffer_[sib].source == trg) {
                // But it is the next element in the buffer.
                ++sib;
                ++trg;
              } else {
                // And we can use it!
                replace = false;
                uint64_t w_i = trg / 64;
                uint64_t w_b = ONE << (trg % 64);
                for (uint16_t bi = 0; bi < 4; ++bi) {
                  if (edges[bi]) {
                    uint64_t* d = new_bits_[bi].data();
#pragma omp atomic
                    d[w_i] = d[w_i] | w_b;
                  }
                }
                // Should be fine to set suffix group start, since there is at
                // most one leader for one suffix group so we will not walk this
                // more than once.
                uint64_t* d = suffix_group_starts_.data();
#pragma omp atomic
                d[w_i] = d[w_i] | w_b;
                break;
              }
            }
          } else {
            // We still need ot update the suffix group leader unless this group
            // becomes empty.
            uint64_t sib = buffer_idx + 1;
            for (uint64_t i = 1; i < 4; ++i) {
              uint64_t trg = buffer_[buffer_idx].source + i;
              if (trg >= n_nodes_) [[unlikely]] {
                // Out of sbwt.
                break;
              }
              if (not suffix_group_starts_[trg]) [[unlikely]] {
                // Next element in sbwt is in same suffix group.
                if (buffer_.size() <= sib || buffer_[sib].source != trg) {
                  // And will not be reomved. So can become new group leader.
                  uint64_t* d = suffix_group_starts_.data();
                  uint64_t w_i = trg / 64;
                  uint64_t w_b = ONE << (trg % 64);
#pragma omp atomic
                  d[w_i] = d[w_i] | w_b;
                } else {
                  // But next element gets removed also.
                  ++sib;
                }
              } else {
                // Next element is in different suffix group.
                break;
              }
            }
          }
          // Replacement was a no go.
          // A new dummy is needed to replace a removed suffix group.
          if (replace) [[unlikely]] {
#pragma omp critical
            { dummies.replace(buffer_[buffer_idx].kmer, edges, *this); }
          }
        }
      }

      if (not buffer_[buffer_idx].b_pred) {
        // Make sure no dangling edge gets left behind.
        uint64_t d_id;
#pragma omp critical
        { d_id = dummies.remove_edge(buffer_[buffer_idx].kmer, *this); }
        // if d_id > 0, parent is a dummy and edge got removed.
        if (d_id == 0) {
          // Parent is not a dummy.
          auto P =
              search_kmer<Kmer_t, k - 1, 0, false>(buffer_[buffer_idx].kmer);
          uint16_t v = buffer_[buffer_idx].kmer.get_v(k - 1);
          uint64_t w_i = P.first / 64;
          uint64_t w_b = ONE << (P.first % 64);
          uint64_t* d = new_bits_[v].data();
#pragma omp atomic
          d[w_i] = d[w_i] | w_b;
        }
      }
    }

    std::vector<uint64_t> removables;
    std::vector<typename decltype(dummies)::Dummy_trie_node> addables;
    dummies.update_removals(removables, addables, *this);

#ifdef VERBOSE
    std::cerr << "Removable dummy indexes: " << std::endl;
    for (auto rd : removables) {
      std::cerr << rd << std::endl;
    }

    std::cerr << "Addable dummies: " << std::endl;
    for (auto ad : addables) {
      std::cerr << ad.to_string() << std::endl;
    }
#endif

#ifdef VERBOSE
    uint64_t xor_count = 0;
    for (uint64_t i = 0; i < n_nodes_; ++i) {
      for (uint16_t ci = 0; ci < 4; ++ci) {
        xor_count += new_bits_[ci][i];
        std::cout << " " << new_bits_[ci][i];
      }
      std::cout << std::endl;
    }
    std::cout << "Total " << xor_count << " xorrable bits" << std::endl;
#endif

#pragma omp parallel for
    for (uint16_t i = 0; i < 4; ++i) {
      new_bits_[i] ^= bits_[i];
    }

#ifdef VERBOSE
    std::cout << "Xorred:\n";
    for (uint64_t i = 0; i < n_nodes_; ++i) {
      for (uint16_t ci = 0; ci < 4; ++ci) {
        std::cout << " " << new_bits_[ci][i];
      }
      std::cout << std::endl;
    }  //*/
#endif

    sdsl::bit_vector old_starts = suffix_group_starts_;

    uint64_t o_size = n_nodes_;
    uint64_t n_size =
        o_size - buffer_.size() + addables.size() - removables.size();
#ifdef DEBUG
    if (n_size > o_size + addables.size() || n_size < 1) {
      std::cerr << "Dummy trie size: " << dummies.nodes.size() << std::endl;
      std::cerr << "Invalid new size\n"
                << "old_size - buffer_size + addable_dummies - "
                << "removabele_dummies = \n"
                << o_size << " - " << buffer_.size() << " + " << addables.size()
                << " - " << removables.size() << " = " << n_size << std::endl;
      exit(1);
    }
#endif
    for (uint16_t i = 0; i < 4; ++i) {
      bits_[i].bit_resize(n_size);
    }
    suffix_group_starts_.bit_resize(n_size);

#ifdef DEBUG
    uint64_t total_compute = 0;
    for (uint16_t ci = 0; ci < 4; ++ci) {
      for (uint64_t i = 0; i < n_nodes_; ++i) {
        total_compute += new_bits_[ci][i];
      }
      for (auto be : buffer_) {
        total_compute -= new_bits_[ci][be.source];
      }
      for (auto de : addables) {
        total_compute += de.out[ci] > 0;
      }
      for (auto de : removables) {
        total_compute -= new_bits_[ci][de];
      }
    }
    if (total_compute != n_size - 1) {
      std::cerr << "Pre-merge edge count missmatch, should be " << n_size - 1
                << " is " << total_compute << std::endl;
    }
#endif

    n_nodes_ = n_size;
    n_kmers_ -= buffer_.size();
    uint64_t s_lim = std::max(n_size, o_size) + 1;

    uint64_t buffer_index = 0;
    uint64_t sbwt_index = 0;
    uint64_t dummy_index = 0;
    uint64_t r_dummy_index = 0;
    uint64_t write_index = 0;
    uint64_t next_i = std::min(
        {buffer_.size() > 0 ? buffer_[buffer_index].source : s_lim,
         addables.size() > 0 ? addables[dummy_index].sbwt_index : s_lim,
         removables.size() > 0 ? removables[r_dummy_index] : s_lim});
    while (sbwt_index < o_size || dummy_index < addables.size()) {
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
              {buffer_.size() > buffer_index ? buffer_[buffer_index].source
                                             : s_lim,
               addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                             : s_lim,
               removables.size() > r_dummy_index ? removables[r_dummy_index]
                                                 : s_lim});
          continue;
        }
      }
      if (r_dummy_index < removables.size() &&
          removables[r_dummy_index] == sbwt_index) {
        ++r_dummy_index;
        ++sbwt_index;
        next_i = std::min(
            {buffer_.size() > buffer_index ? buffer_[buffer_index].source
                                           : s_lim,
             addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                           : s_lim,
             removables.size() > r_dummy_index ? removables[r_dummy_index]
                                               : s_lim});
        continue;
      }
      if (buffer_index < buffer_.size() &&
          buffer_[buffer_index].source == sbwt_index) {
        ++buffer_index;
        ++sbwt_index;
        next_i = std::min(
            {buffer_.size() > buffer_index ? buffer_[buffer_index].source
                                           : s_lim,
             addables.size() > dummy_index ? addables[dummy_index].sbwt_index
                                           : s_lim,
             removables.size() > r_dummy_index ? removables[r_dummy_index]
                                               : s_lim});
        continue;
      }
      uint64_t copy_count = std::min(next_i, o_size) - sbwt_index;
      copy_count = std::min(n_size - write_index, copy_count);
#ifdef DEBUG
      if (sbwt_index + copy_count > o_size) {
        std::cerr << "Can't copy " << copy_count << " elements from "
                  << sbwt_index << "\nGoes " << sbwt_index + copy_count - o_size
                  << std::endl;
        exit(1);
      }
      if (write_index + copy_count > n_size) {
        std::cerr << "Can't copy " << copy_count << " elements to "
                  << write_index << "\nGoes "
                  << write_index + copy_count - n_size << std::endl;
        exit(1);
      }
#endif
      uint64_t w;
      uint64_t w_i = write_index / 64;
      uint8_t w_o = write_index % 64;
      write_index += copy_count;
      uint64_t r_i = sbwt_index / 64;
      uint8_t r_o = sbwt_index % 64;
      std::array<uint64_t*, 4> w_ptr = {
          bits_[0].data() + w_i, bits_[1].data() + w_i, bits_[2].data() + w_i,
          bits_[3].data() + w_i};
      uint64_t* s_w_ptr = suffix_group_starts_.data() + w_i;
      std::array<uint64_t*, 4> r_ptr = {
          new_bits_[0].data() + r_i, new_bits_[1].data() + r_i,
          new_bits_[2].data() + r_i, new_bits_[3].data() + r_i};
      uint64_t* s_r_ptr = old_starts.data() + r_i;
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
      sbwt_index = next_i;
    }

    for (uint16_t i = 0; i < 4; ++i) {
      sdsl::util::init_support(rank_supports_[i], &(bits_[i]));
    }

    for (uint16_t i = 1; i < 4; ++i) {
      c_table_[i] = rank_supports_[i - 1].rank(n_size);
      c_table_[i] += c_table_[i - 1];
    }

#ifdef VERBOSE
    print();
#endif

    compute_k_mer_precalc();
    buffer_.clear();
    for (uint16_t i = 0; i < 4; ++i) {
      new_bits_[i].bit_resize(0);
    }
#ifdef DEBUG
    if (not is_valid()) {
      std::cerr << "Validation fail!" << std::endl;
    }
#endif
  }

 public:
  Buffered_SBWT(double buffer_gigs = 0.5)
      : bits_(),
        rank_supports_(),
        new_bits_(),
        suffix_group_starts_(1),
        c_table_({1, 1, 1, 1}),
        kmer_prefix_precalc_(ONE << (2 * precalc_k), {1, 1}),
        n_nodes_(1),
        n_kmers_(0) {
    for (uint16_t i = 0; i < 4; ++i) {
      bits_[i].resize(1);
      sdsl::util::init_support(rank_supports_[i], &bits_[i]);
    }
    suffix_group_starts_[0] = true;
    buffer_limit_ = gigs_to_elems<B_elem>(buffer_gigs);
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

  template <class OS>
  uint64_t serialize(OS& out_stream) const {
    if constexpr (std::is_same_v<std::string, OS>) {
      throwing_ofstream out(out_stream, ios::binary);
      return serialize(out);
    } else {
      uint64_t written = 0;
      written += serialize_string(SBWT_VERSION, out_stream);
      for (uint16_t i = 0; i < 4; ++i) {
        written += bits_[i].serialize(out_stream.stream);
        written += rank_supports_[i].serialize(out_stream.stream);
      }
      written += suffix_group_starts_.serialize(out_stream.stream);
      written += serialize_std_array(c_table_, out_stream);
      written += serialize_std_vector(kmer_prefix_precalc_, out_stream);
      out_stream.write(reinterpret_cast<const char*>(&n_nodes_),
                       sizeof(n_nodes_));
      written += sizeof(n_nodes_);
      out_stream.write(reinterpret_cast<const char*>(&n_kmers_),
                       sizeof(n_kmers_));
      written += sizeof(n_kmers_);
      uint16_t l_k = k;
      out_stream.write(reinterpret_cast<char*>(&l_k), sizeof(l_k));
      written += sizeof(l_k);
      return written;
    }
  }

  template <class OS>
  uint64_t serialize_old_format(OS& out_stream) const {
    if constexpr (std::is_same_v<std::string, OS>) {
      throwing_ofstream out(out_stream, ios::binary);
      return serialize_old_format(out);
    } else {
      uint64_t written = 0;
      written += serialize_string(SBWT_VARIANT, out_stream);
      std::string old_version = "v0.1";
      written += serialize_string(old_version, out_stream);
      for (uint16_t i = 0; i < 4; ++i) {
        bits_[i].serialize(out_stream.stream);
      }
      for (uint16_t i = 0; i < 4; ++i) {
        rank_supports_[i].serialize(out_stream.stream);
      }
      written += suffix_group_starts_.serialize(out_stream.stream);
      written += serialize_std_array(c_table_, out_stream);
      std::vector<std::pair<int64_t, int64_t>> old_precalc;
      for (auto P : kmer_prefix_precalc_) {
        std::pair<int64_t, int64_t> upd_P(P);
        std::pair<int64_t, int64_t> unfound_P = {-1, -1};
        old_precalc.push_back(P.first == P.second ? unfound_P : upd_P);
      }
      written += serialize_std_vector(old_precalc, out_stream);
      int64_t writable_v = precalc_k;
      out_stream.write(reinterpret_cast<char*>(&writable_v),
                       sizeof(writable_v));
      written += sizeof(writable_v);

      writable_v = n_nodes_;
      out_stream.write(reinterpret_cast<char*>(&writable_v),
                       sizeof(writable_v));
      written += sizeof(writable_v);

      writable_v = n_kmers_;
      out_stream.write(reinterpret_cast<char*>(&writable_v),
                       sizeof(writable_v));
      written += sizeof(writable_v);

      writable_v = k;
      out_stream.write(reinterpret_cast<char*>(&writable_v),
                       sizeof(writable_v));
      written += sizeof(writable_v);

      return written;
    }
  }

  template <class IS>
  void load(IS& in_stream) {
    if constexpr (std::is_same_v<std::string, IS>) {
      throwing_ifstream in(in_stream, ios::binary);
      load(in);
    } else {
      std::string var = load_string(in_stream);
      if (var == SBWT_VERSION) {
        for (uint16_t i = 0; i < 4; ++i) {
          bits_[i].load(in_stream.stream);
          rank_supports_[i].load(in_stream.stream, &(bits_[i]));
        }
        suffix_group_starts_.load(in_stream.stream);
        load_std_array(c_table_, in_stream);
        load_std_vector(kmer_prefix_precalc_, in_stream);
        in_stream.read(reinterpret_cast<char*>(&n_nodes_), sizeof(n_nodes_));
        in_stream.read(reinterpret_cast<char*>(&n_kmers_), sizeof(n_kmers_));
        uint16_t t_k;
        in_stream.read(reinterpret_cast<char*>(&t_k), sizeof(t_k));
        if (t_k != k) {
          throw std::runtime_error(
              "Attempt to load sbwt with the wrong value for k");
        }
      } else if (var == SBWT_VARIANT) {
        var = load_string(in_stream);
        if (var == "v0.1") {
          for (uint16_t i = 0; i < 4; ++i) {
            bits_[i].load(in_stream.stream);
          }
          for (uint16_t i = 0; i < 4; ++i) {
            rank_supports_[i].load(in_stream.stream, &(bits_[i]));
          }
          suffix_group_starts_.load(in_stream.stream);
          load_std_array(c_table_, in_stream);
          load_std_vector(kmer_prefix_precalc_, in_stream);
          int64_t t_p_k;
          in_stream.read(reinterpret_cast<char*>(&t_p_k), sizeof(t_p_k));
          in_stream.read(reinterpret_cast<char*>(&n_nodes_), sizeof(n_nodes_));
          in_stream.read(reinterpret_cast<char*>(&n_kmers_), sizeof(n_kmers_));
          int64_t t_k;
          in_stream.read(reinterpret_cast<char*>(&t_k), sizeof(t_k));
          if (t_k != k) {
            throw std::runtime_error(
                "Attempt to load sbwt with the wrong value for k");
          }
          if (t_p_k != precalc_k) {
            throw std::runtime_error(
                "Attempt to load sbwt with the wrong value for precalc_k");
          }
          compute_k_mer_precalc();
        } else {
          throw std::runtime_error("Invalid version string or corrupted SBWT");
        }
      } else {
        throw std::runtime_error("Invalid variant/version or corrupted SBWT");
      }
    }
  }

  void add(const std::string& addable) {
    add_string(addable, buffer_);
    if (buffer_.size() == 0) {
      return;
    }
    commit();
  }

  void del(const std::string& removable) {
    delete_string(removable, buffer_);
    if (buffer_.size() == 0) {
      return;
    }
    r_commit();
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
                buffer_.insert(buffer_.end(), addable_k_mers.begin(),
                               addable_k_mers.end());
              }
            }
            break;
          }
          add_string(addable, addable_k_mers);
#pragma omp atomic
          k_mer_count = k_mer_count + addable_k_mers.size();

          if (k_mer_count >= buffer_limit_) {
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
                buffer_.insert(buffer_.end(), addable_k_mers.begin(),
                               addable_k_mers.end());
              }
            }
            break;
          }
        }
      }
      if (buffer_.size() >= buffer_limit_) {
        commit();
      }
    }
    if (buffer_.size() > 0) {
      commit();
    }
    return total_k_mers;
  }

  template <class removables_t>
  uint64_t del_all(removables_t& addables) {
    uint64_t total_k_mers = 0;
    bool cont = true;
    while (cont) {
      uint64_t k_mer_count = 0;
// parallel
#pragma omp parallel
      {
        std::vector<B_elem> removable_k_mers;
        while (true) {
          bool do_break = false;
          std::string removable;
#pragma omp critical
          {
            do_break = not addables.get(removable);
            cont = cont && not do_break;
            if (not do_break) {
              total_k_mers += removable.size() - k + 1;
            }
          }
          if (do_break) {
            if (removable_k_mers.size() > 0) {
#pragma omp critical
              {
                buffer_.insert(buffer_.end(), removable_k_mers.begin(),
                               removable_k_mers.end());
              }
            }
            break;
          }
          delete_string(removable, removable_k_mers);
#pragma omp atomic
          k_mer_count = k_mer_count + removable_k_mers.size();

          if (k_mer_count >= buffer_limit_) {
            if (removable_k_mers.size() > 0) {
              std::sort(removable_k_mers.begin(), removable_k_mers.end());
              uint64_t trg = 0;
              for (uint64_t i = 0; i < removable_k_mers.size(); ++i) {
                if (removable_k_mers[trg].kmer != removable_k_mers[i].kmer) {
                  removable_k_mers[++trg] = removable_k_mers[i];
                }
              }
              ++trg;
              removable_k_mers.resize(trg);
#pragma omp critical
              {
                buffer_.insert(buffer_.end(), removable_k_mers.begin(),
                               removable_k_mers.end());
              }
            }
            break;
          }
        }
      }
      if (buffer_.size() >= buffer_limit_) {
        r_commit();
      }
    }
    if (buffer_.size() > 0) {
      r_commit();
    }
    return total_k_mers;
  }

  uint64_t number_of_kmers() const { return n_kmers_; }

  int64_t search(const std::string& kmer) const {
    Kmer_t km(kmer);
    auto I = search_kmer(km);
    return I.first == I.second ? -1 : I.first;
  }

  void streaming_search(const std::string& input, uint64_t& count,
                        uint64_t& matches) {
    Streaming_search<true> ss(input, *this);
    uint64_t i;
    bool found;
    while (ss.next(i, found)) {
      matches += found;
      ++count;
    }
  }

  std::vector<bool> streaming_search(const std::string& input) const {
    std::vector<bool> ret;
    Streaming_search<true> ss(input, *this);
    uint64_t i;
    bool found;
    while (ss.next(i, found)) {
      ret.push_back(found);
    }
    return ret;
  }

  uint16_t get_precalc_k() const { return precalc_k; }

  uint16_t get_k() const { return k; }

  std::vector<std::string> reconstruct_all_kmers() const {
    uint64_t n_nodes = n_nodes_;
    std::vector<uint64_t> C_array(4);

    std::vector<char> last;  // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for (uint64_t i = 0; i < n_nodes; i++)
      if (bits_[0][i]) last.push_back('A');

    C_array[1] = last.size();
    for (uint64_t i = 0; i < n_nodes; i++)
      if (bits_[1][i]) last.push_back('C');

    C_array[2] = last.size();
    for (uint64_t i = 0; i < n_nodes; i++)
      if (bits_[2][i]) last.push_back('G');

    C_array[3] = last.size();
    for (uint64_t i = 0; i < n_nodes; i++)
      if (bits_[3][i]) last.push_back('T');

    if (last.size() != n_nodes) {
      cerr << "BUG " << last.size() << " " << n_nodes << endl;
      exit(1);
    }

    std::string kmers_concat(n_nodes * k, '\0');

    for (uint64_t round = 0; round < k; round++) {
      // cerr << "round " << round << "/" << k-1 << endl;
      for (uint64_t i = 0; i < n_nodes; i++) {
        uint64_t pos = k - 1 - round;
        kmers_concat[i * k + pos] = last[i];
      }

      // Propagate the labels one step forward in the graph
      std::vector<char> propagated(n_nodes, '$');
      uint64_t A_ptr = C_array[0];
      uint64_t C_ptr = C_array[1];
      uint64_t G_ptr = C_array[2];
      uint64_t T_ptr = C_array[3];
      for (uint64_t i = 0; i < n_nodes; i++) {
        if (bits_[0][i]) propagated[A_ptr++] = last[i];
        if (bits_[1][i]) propagated[C_ptr++] = last[i];
        if (bits_[2][i]) propagated[G_ptr++] = last[i];
        if (bits_[3][i]) propagated[T_ptr++] = last[i];
      }
      last = propagated;
    }
    std::vector<std::string> ret;
    for (uint64_t i = 0; i < kmers_concat.size() - k + 1; i += k) {
      ret.push_back(kmers_concat.substr(i, k));
    }
    return ret;
  }

  void print() {
    std::cout << "SBWT C-table, matrix, suffix group starts." << std::endl;
    for (uint16_t i = 0; i < 4; ++i) {
      std::cout << c_table_[i] << " ";
    }
    std::cout << "<- C-table" << std::endl;
    for (uint16_t i = 0; i <= k; ++i) {
      std::cout << " ";
    }
    std::cout << "a c g t    sgs" << std::endl;
    auto r = reconstruct_all_kmers();
    for (uint64_t i = 0; i < n_nodes_; ++i) {
      std::cout << r[i] << " " << bits_[0][i] << " " << bits_[1][i] << " "
                << bits_[2][i] << " " << bits_[3][i] << "    "
                << suffix_group_starts_[i] << std::endl;
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
    c = 0;
    for (uint16_t i = 0; i < 4; ++i) {
      uint64_t rank_count = rank_supports_[i].rank(n);
      uint64_t raw_count = 0;
      for (uint64_t ii = 0; ii < n; ++ii) {
        raw_count += bits_[i][ii];
      }
      if (rank_count != raw_count) {
        std::cerr << i << " bv rank missmatch: rank(" << n
                  << ") = " << rank_count << ", should be " << raw_count
                  << std::endl;
        ok = false;
      }
      c += rank_count;
    }
    if (c != n - 1) {
      std::cerr << "Tree is not. number of edges should be " << n - 1 << ", is "
                << c << std::endl;
      ok = false;
    }
    return ok;
  }

  bool compare(const Buffered_SBWT& other) const {
    bool is_ok = true;
    if (n_kmers_ != other.n_kmers_) {
      std::cerr << "k-mer count missmatch: " << n_kmers_
                << " != " << other.n_kmers_ << std::endl;
      is_ok = false;
    }
    if (n_nodes_ != other.n_nodes_) {
      std::cerr << "node count (dummy count) missmatch: " << n_nodes_
                << " != " << other.n_nodes_ << std::endl;
      is_ok = false;
    }
    const std::array<char, 4> alpha = {'A', 'C', 'G', 'T'};
    for (uint16_t i = 0; i < 4; ++i) {
      if (bits_[i].size() != other.bits_[i].size()) {
        std::cerr << alpha[i] << " bv size missmatch\n    " << bits_[i].size()
                  << " != " << other.bits_[i].size() << std::endl;
        is_ok = false;
      }
      for (uint64_t bv_i = 0;
           bv_i < std::min(bits_[i].size(), other.bits_[i].size()); ++bv_i) {
        if (bits_[i][bv_i] != other.bits_[i][bv_i]) {
          std::cerr << alpha[i] << " bv content missmatch at " << bv_i
                    << std::endl;
          is_ok = false;
          break;
        }
      }
    }
    if (suffix_group_starts_.size() != other.suffix_group_starts_.size()) {
      std::cerr << "Suffix group start bv size missmatch\n    "
                << suffix_group_starts_.size()
                << " != " << other.suffix_group_starts_.size() << std::endl;
      is_ok = false;
    }
    if (suffix_group_starts_ != other.suffix_group_starts_) {
      std::cerr << "Suffix group start bv content missmatch" << std::endl;
      is_ok = false;
    }
    if (kmer_prefix_precalc_.size() != other.kmer_prefix_precalc_.size()) {
      std::cerr << "Precalc kmer count missmatch\n   "
                << kmer_prefix_precalc_.size()
                << " != " << other.kmer_prefix_precalc_.size() << std::endl;
    }
    uint16_t r_count = 0;
    uint16_t e_limit = 10;
    for (uint64_t i = 0; i < kmer_prefix_precalc_.size(); ++i) {
      if (r_count < e_limit &&
          kmer_prefix_precalc_[i] != other.kmer_prefix_precalc_[i]) {
        std::cerr << "Precalc kmer " << i << " missmatch:\n"
                  << "(" << kmer_prefix_precalc_[i].first << ", "
                  << kmer_prefix_precalc_[i].second << ") <-> ("
                  << other.kmer_prefix_precalc_[i].first << ", "
                  << other.kmer_prefix_precalc_[i].second << ")" << std::endl;
        is_ok = false;
        ++r_count;
        if (r_count == e_limit) {
          std::cerr << "Reported " << e_limit
                    << " errors before not looking anymore" << std::endl;
        }
      }
    }
    return is_ok;
  }
};

template <uint16_t k, uint16_t precalc_k>
template <bool short_circuit>
class Buffered_SBWT<k, precalc_k>::Streaming_search {
 private:
  const Buffered_SBWT& static_bsbwt_;
  const std::string& searchable_;
  Kmer<k> kmer_;
  uint64_t pred_;
  uint64_t loc_ = 0;
  uint64_t index_ = k - 1;
  uint16_t p_size_ = 0;
  bool found_ = false;

 public:
  Streaming_search(const std::string& searchable, const Buffered_SBWT& bsbwt)
      : static_bsbwt_(bsbwt), searchable_(searchable), kmer_(searchable) {}

  bool next(uint64_t& idx, bool& found) {
    if (index_ >= searchable_.size()) {
      return false;
    }
    if (static_bsbwt_.n_nodes_ == 1) [[unlikely]] {
      loc_ = 1;
      found_ = false;
      pred_ = 1;
      p_size_ = 0;
    } else if (found_) {
      uint64_t a = loc_;
      uint64_t b = a + 1;
      for (uint16_t i = 0; i < 4; ++i) {
        if (static_bsbwt_.suffix_group_starts_[a]) {
          break;
        }
        --a;
      }
      pred_ = a;
      p_size_ = b - a;
      static_bsbwt_.fl(a, b, kmer_.get_v(k - 1));
      loc_ = a;
      found_ = a < b;
    } else {
      auto precalc_idx = kmer_.template get_first_v<precalc_k>();
      auto p = static_bsbwt_.kmer_prefix_precalc_[precalc_idx];
      uint64_t a = p.first;
      uint64_t b = p.second;
      for (uint16_t i = precalc_k; i < k - 1; ++i) {
        static_bsbwt_.fl(a, b, kmer_.get_v(i));
        if constexpr (short_circuit) {
          if (a >= b) {
            found_ = false;
            found = found_;
            ++index_;
            if (index_ < searchable_.size()) [[likely]] {
              kmer_ = kmer_.next(searchable_[index_]);
            }
            return true;
          }
        }
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
            Kmer<k>& km) {
    km = kmer_;
    bool r = next(idx, found);
    pred = pred_;
    pred_size = p_size_;
    return r;
  }

  bool next(uint64_t& idx, bool& found, Kmer<k>& km) {
    km = kmer_;
    bool r = next(idx, found);
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

    std::string to_string() {
      std::string ret = std::to_string(kmer_length);
      ret.append(":").append(kmer.to_string());
      for (auto o : out) {
        ret.append(" ").append(std::to_string(o));
      }
      ret.append(", ").append(std::to_string(sbwt_index));
      return ret;
    }
  };

  std::vector<Dummy_trie_node> nodes;

 private:
  void walk(uint64_t index, uint16_t depth, const Buffered_SBWT& sbwt) {
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
  Dummy_trie(const Buffered_SBWT& sbwt) : nodes() {
    nodes.push_back({Kmer<k>(), {0, 0, 0, 0}, 0, 0, true});
    walk(0, 1, sbwt);
  }

  template <class Kmer_t>
  uint64_t mark_for_removal(const Kmer_t& kmer) {
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
  uint64_t remove_edge(const Kmer_t& kmer, Buffered_SBWT& sbwt) {
    uint64_t dummy_id = 0;
    for (uint16_t i = 0; i < k - 1; ++i) {
      uint16_t v = kmer.get_v(i);
      dummy_id = nodes[dummy_id].out[v];
      if (dummy_id == 0) {
        return 0;
      }
    }
    uint16_t v = kmer.get_v(k - 1);
    nodes[dummy_id].out[v] = 0;
    if (nodes[dummy_id].sbwt_index > 0) {
      sbwt.new_bits_[v][nodes[dummy_id].sbwt_index] = true;
    }
    nodes[dummy_id].keep = nodes[dummy_id].out[0] > 0;
    for (uint16_t i = 1; i < 4; ++i) {
      nodes[dummy_id].keep |= nodes[dummy_id].out[i] > 0;
    }
    return dummy_id;
  }

  template <class Kmer_t, class arr_t>
  void replace(const Kmer_t& kmer, const arr_t& edges, Buffered_SBWT& sbwt) {
    uint64_t id = 0;
    uint64_t depth = 1;
    while (depth < k) {
      uint16_t v = kmer.get_v(depth);
      if (nodes[id].out[v] > 0) {
        id = nodes[id].out[v];
      } else {
        if (id == 0 || nodes[id].sbwt_index != 0) {
          sbwt.new_bits_[v][nodes[id].sbwt_index] = true;
        }
        break;
      }
      ++depth;
    }
    while (depth < k) {
      uint64_t idx = nodes.size();
      uint16_t v = kmer.get_v(depth);
      Dummy_trie_node nd = {
          nodes[id].kmer.next_v(v), {0, 0, 0, 0}, 0, uint16_t(depth), true};
      if (depth == k - 1) {
        for (uint16_t i = 0; i < 4; ++i) {
          nd.out[i] = edges[i];
        }
      }
      nodes.push_back(nd);
      nodes[id].out[v] = idx;
      id = idx;
      ++depth;
    }
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