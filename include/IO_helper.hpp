#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

#include "SeqIO/SeqIO.hh"

namespace sbwt {

template <uint16_t k, bool gzipped = false>
class io_container {
 private:
  typedef std::conditional<
      gzipped,
      seq_io::Multi_File_Reader<seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>>>,
      seq_io::Multi_File_Reader<>>::type reader_t;

  reader_t reader;
  bool filter_n_;
  uint64_t s_len;
  uint64_t s_loc;

 public:
  io_container(std::vector<std::string>& filelist, bool rev_comp, bool filter_n)
      : reader(filelist), filter_n_(filter_n) {
    if (rev_comp) reader.enable_reverse_complements();
    s_len = reader.get_next_read_to_buffer();
    s_loc = 0;
  }

  bool get(std::string& out) {
    while (true) {
      if (s_len == 0) {
        out.clear();
        return false;
      }
      if (s_len <= s_loc) {
        s_len = reader.get_next_read_to_buffer();
        s_loc = 0;
        if (s_len == 0) {
          out.clear();
          return false;
        }
      }
      uint64_t end = s_len;
      if (filter_n_) {
        while (s_loc < s_len) {
          char c = reader.read_buf[s_loc];
          if (c == 'A' || c == 'C' || c == 'G' || c == 'T') [[likely]] {
            break;
          }
          ++s_loc;
        }
        end = s_loc;
        while (end < s_len) {
          char c = reader.read_buf[end];
          if (c == 'A' || c == 'C' || c == 'G' || c == 'T') [[likely]] {
            ++end;
          } else {
            break;
          }
        }
      }
      if (end - s_loc >= k) {
        out = std::string(reader.read_buf + s_loc, end - s_loc);
        s_loc = end;
        return true;
      }
      s_loc = end;
    }
  }
};
}  // namespace sbwt