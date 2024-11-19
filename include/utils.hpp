#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>

namespace sbwt {

template <class out_t>
inline uint64_t serialize_string(const std::string& S, out_t& out) {
  uint64_t size = S.size();
  out.write(reinterpret_cast<char*>(&size), sizeof(size));
  out.write(S.data(), size);
  return sizeof(size) + size;
}

template <class in_t>
inline std::string load_string(in_t& in) {
  uint64_t size;
  in.read(reinterpret_cast<char*>(&size), sizeof(size));
  string S(size, '\0');
  in.read(S.data(), size);
  return S;
}

// Utility function: Serialization for a std::vector
// Returns number of bytes written
template <typename T, class out_t>
inline uint64_t serialize_std_vector(const std::vector<T>& v, out_t& os) {
  uint64_t n_bytes = sizeof(T) * v.size();
  os.write(reinterpret_cast<char*>(&n_bytes), sizeof(n_bytes));
  os.write(reinterpret_cast<char*>(v.data()), n_bytes);
  return sizeof(n_bytes) + n_bytes;
}

template <typename T, class in_t>
inline void load_std_vector(std::vector<T>& v, in_t& is) {
  uint64_t n_bytes = 0;
  is.read(reinterpret_cast<char*>(&n_bytes), sizeof(n_bytes));
  assert(n_bytes % sizeof(T) == 0);
  if (n_bytes / sizeof(T) != v.size()) {
    throw std::runtime_error("Invalid vector size");
  }
  is.read(reinterpret_cast<char*>(v.data()), n_bytes);
}

// Utility function: Serialization for a std::vector
// Returns number of bytes written
template <typename AT, class out_t>
inline uint64_t serialize_std_array(const AT& a, out_t& os) {
  uint64_t n_bytes = sizeof(AT);
  os.write(reinterpret_cast<char*>(&n_bytes), sizeof(n_bytes));
  os.write(reinterpret_cast<char*>(a.data()), sizeof(AT));
  return sizeof(AT) + sizeof(n_bytes);
}

template <typename AT, class in_t>
inline void load_std_array(AT& a, in_t& is) {
  uint64_t n_bytes = 0;
  is.read(reinterpret_cast<char*>&n_bytes, sizeof(n_bytes));
  assert(n_bytes % sizeof(T) == 0);
  if (n_bytes / sizeof(T) != a.size()) {
    throw std::runtime_error("Invalid array size");
  }
  is.read(reinterpret_cast<char*>(a.data()), sizeof(AT));
}

}  // namespace sbwt