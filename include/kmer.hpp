#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace sbwt {

const static constexpr std::array<char, 4> from_bits_to_char_table_ = {
    'A', 'C', 'G', 'T'};

const static constexpr std::array<uint8_t, 256> from_char_to_bits_table = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

template <uint16_t len>
class kmer_colex_compare;  // Compiler needs this forward declaration to make it
                           // a friend class

// A bit-packed k-mer class containing at most max_len characters
// Maximum length supported is 255 because the length is stored in an uint8_t
template <uint16_t len>
class Kmer {
 public:
  friend class kmer_colex_compare<len>;
  const static constexpr uint16_t KMER_LEN = len;
 private:
  // Two bits per character.
  // The rightmost character of the whole k-mer is at the most significant
  // bits of block 0. See the function `get` for a precise definition of the
  // data layout.
  std::array<uint64_t, (len + 31) / 32> data_;

  char constexpr to_char(uint8_t x) const {
    assert(x < 4);
    return from_bits_to_char_table_[x];
  }

  uint8_t constexpr to_bitpair(char c) const {
    return from_char_to_bits_table[c];
  }

 public:
  Kmer() : data_() {}

  Kmer(const std::string& S) : data_() {
    for (uint64_t i = 0; i < len; i++) {
      set(i, S[i]);
    }
  }

  Kmer(const char* S) : data_() {
    for (int64_t i = 0; i < len; i++) {
      set(i, S[i]);
    }
  }

  Kmer next(char c) { return next_v(to_bitpair(c)); }

  Kmer next_v(uint64_t v) {
    Kmer ret(*this);
    for (uint16_t i = ret.data_.size() - 1; i > 0; --i) {
      ret.data_[i] >>= 2;
      ret.data_[i] |= ret.data_[i - 1] << 62;
    }
    ret.data_[0] >>= 2;
    ret.data_[0] |= v << 62;
    if constexpr (len % 32) {
      const constexpr uint16_t lwb = (len % 32) * 2;
      const constexpr uint64_t mask = ~((uint64_t(1) << (64 - lwb)) - 1);
      ret.data_[ret.data_.size() - 1] &= mask;
    }
    return ret;
  }

  Kmer& operator++() {
    const constexpr uint64_t v = uint64_t(1) << (2 * ((32 - (len % 32)) % 32));
    for (uint16_t i = 0; i < data_.size() - 1; ++i) {
      ++data_[i];
      if (data_[i] != 0) [[likely]] {
        return *this;
      }
    }
    data_[data_.size() - 1] += v;
    return *this;
  }

  void clear() { data_.clear(); }

  // Get the character at index idx from the start
  char get(int64_t idx) const { return to_char(get_v(idx)); }

  uint16_t get_v(uint64_t idx) const {
    assert(idx < len);
    uint64_t block_idx = 0;
    if constexpr (len > 32) {
      block_idx = (len - 1 - idx) / 32;
    }
    uint64_t block_offset = (len - 1 - idx) % 32;
    return (data_[block_idx] >> ((31 - block_offset) * 2)) & 0x03;
  }

  uint64_t get_first_v(uint16_t i) const {
    uint64_t ret = 0;
    for (uint64_t c = 0; c < i; ++c) {
      ret = (ret << 2) | get_v(i - 1 - c);
    }
    return ret;
  }

  // Get the character at index idx from the start
  void set(int64_t idx, char c) {
    assert(idx >= 0 && idx < len);
    assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
    int64_t block_idx = 0;
    if constexpr (len > 32) {
      block_idx = (len - 1 - idx) / 32;
    }

    int64_t block_offset = (len - 1 - idx) % 32;
    data_[block_idx] &= (~(uint64_t)0) ^
                        ((uint64_t)0x03 << ((31 - block_offset) * 2));  // Clear
    data_[block_idx] |= (uint64_t)to_bitpair(c)
                        << ((31 - block_offset) * 2);  // Set
  }

  bool operator==(const Kmer& other) const {
    bool ret = this->data_[0] == other.data_[0];
    for (int64_t i = 1; i < data_.size(); i++) {
      ret &= this->data_[i] == other.data_[i];
    }
    return ret;
  }

  bool operator!=(const Kmer& other) const { return !(*this == other); }

  // Strict colexicographic comparison
  // Returns true if this is smaller than the other
  bool operator<(const Kmer& other) const {
    for (int64_t i = 0; i < data_.size(); i++) {
      if (this->data_[i] < other.data_[i]) return true;
      if (this->data_[i] > other.data_[i]) return false;
    }
    return false;
  }

  bool operator>(const Kmer& other) const {
    return !(*this < other) && !(*this == other);
  }

  bool operator<=(const Kmer& rhs) const { return !(*this > rhs); }

  bool operator>=(const Kmer& rhs) const { return !(*this < rhs); }

  Kmer copy() const {
    Kmer other(*this);
    return other;
  }

  std::string to_string(uint16_t a = 0, uint16_t b = 0) const {
    if (a > b) {
      b = a;
      a = 0;
    } else if (a == 0 && b == 0) {
      b = len;
    }
    std::string S;
    for (int64_t i = a; i < b; i++) S += get(i);
    return S;
  }

  char last() const { return get(len - 1); }

  char first() const { return get(0); }

  static constexpr int64_t size_in_bytes() { return sizeof(Kmer); }

  template<class out_t>
  void serialize(out_t& out) const {
    out.write((char*)(data_.data()), sizeof(data_));
  }

  template<class in_t>
  void load(in_t& in) { in.read((char*)(data_.data()), sizeof(data_)); }

  // `out` must have at least size_in_bytes() bytes of space
  void serialize(char* out) const {
    std::copy(reinterpret_cast<char*>(data_.data()),
              reinterpret_cast<char*>(data_.data() + data_.size()), out);
  }

  // Load from memory serialized to by member function `serialize(char* out)`
  void load(const char* in) {
    std::copy(reinterpret_cast<const uint64_t*>(in),
              reinterpret_cast<const uint64_t*>(in) + data_.size(),
              data_.data());
  }

  static constexpr uint64_t size() { return len; }

  Kmer suf() const {
    auto o = copy();
    if constexpr (len % 32 == 1) {
      o.data_[data_.size() - 1] = 0;
    } else {
      const constexpr uint16_t lwb = 2 * ((len - 1) % 32);
      const constexpr uint64_t mask = ~((uint64_t(1) << (64 - lwb)) - 1);
      o.data_[data_.size() - 1] &= mask;
    }
    return o;
  }

  Kmer pref() const {
    auto o = copy();
    for (uint16_t i = 1; i < data_.size(); ++i) {
      o.data_[i - 1] <<= 2;
      o.data_[i - 1] |= o.data_[i] >> 62;
    }
    o.data_[data_.size() - 1] <<= 2;
    return o;
  }

  uint64_t miw() const { return data_[0]; }
};

// Strict colexicographic comparison
template <uint16_t max_len>
struct kmer_colex_compare {
  inline bool operator()(const Kmer<max_len>& A, const Kmer<max_len>& B) {
    return A < B;
  }
};

}  // namespace sbwt
