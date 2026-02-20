#include "obj_parser.h"

#include <charconv>
#include <fstream>
#include <algorithm>
#include <cctype>

namespace s21 {

namespace {
inline bool IsSpace(char c) {
  return std::isspace(static_cast<unsigned char>(c)) != 0;
}
} // namespace

bool ObjParser::Parse(const std::string& filename,
                      std::vector<Vertex>& vertices,
                      std::vector<Edge>& edges,
                      std::vector<std::uint32_t>& tri_indices) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) return false;

  vertices.clear();
  edges.clear();
  tri_indices.clear();

  // Оценочное резервирование
  file.seekg(0, std::ios::end);
  std::streampos file_size = file.tellg();
  file.seekg(0, std::ios::beg);

  vertices.reserve(static_cast<size_t>(file_size / 30));
  edges.reserve(static_cast<size_t>(file_size / 20));
  tri_indices.reserve(static_cast<size_t>(file_size / 20)); // грубо: 1 триангл ~ 20 байт текста

  std::string line;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;

    // ----- Vertex -----
    if (line[0] == 'v' && line.size() > 1 && IsSpace(line[1])) {
      const char* start = line.data();
      const char* end   = start + line.size();
      const char* ptr   = start + 1;

      while (ptr < end && IsSpace(*ptr)) ++ptr;

      double x = 0.0, y = 0.0, z = 0.0;

      auto [p1, ec1] = std::from_chars(ptr, end, x);
      if (ec1 != std::errc()) continue;
      while (p1 < end && IsSpace(*p1)) ++p1;

      auto [p2, ec2] = std::from_chars(p1, end, y);
      if (ec2 != std::errc()) continue;
      while (p2 < end && IsSpace(*p2)) ++p2;

      auto [p3, ec3] = std::from_chars(p2, end, z);
      if (ec3 != std::errc()) continue;

      vertices.push_back({x, y, z});
      continue;
    }

    // ----- Face -----
    if (line[0] == 'f' && line.size() > 1 && IsSpace(line[1])) {
      const char* start = line.data();
      const char* end   = start + line.size();
      const char* ptr   = start + 1;

      while (ptr < end && IsSpace(*ptr)) ++ptr;

      std::vector<int> indices;
      indices.reserve(16);

      bool bad_face = false;

      while (ptr < end) {
        while (ptr < end && IsSpace(*ptr)) ++ptr;
        if (ptr >= end) break;

        // комментарий в конце строки
        if (*ptr == '#') break;

        const char* token_start = ptr;
        while (ptr < end && !IsSpace(*ptr)) ++ptr;
        const char* token_end = ptr;

        const char* idx_end = token_start;
        while (idx_end < token_end && *idx_end != '/') ++idx_end;

        if (idx_end == token_start) { bad_face = true; break; }

        int idx = 0;
        auto [p, ec] = std::from_chars(token_start, idx_end, idx);
        if (ec != std::errc()) { bad_face = true; break; }

        // 1-base / negative -> 0-base
        if (idx < 0) {
          idx = static_cast<int>(vertices.size()) + idx; // idx отрицательный
        } else {
          idx -= 1;
        }

        // Валидация диапазона (критично для безопасности)
        if (idx < 0 || idx >= static_cast<int>(vertices.size())) {
          bad_face = true;
          break;
        }

        indices.push_back(idx);
      }

      if (bad_face || indices.size() < 3) continue;

      for (size_t i = 0; i < indices.size(); ++i) {
        size_t j = (i + 1) % indices.size();
        unsigned int v1 = static_cast<unsigned int>(indices[i]);
        unsigned int v2 = static_cast<unsigned int>(indices[j]);
        if (v1 == v2) continue;
        if (v1 > v2) std::swap(v1, v2);
        edges.emplace_back(v1, v2);
      }

      const std::uint32_t v0 = static_cast<std::uint32_t>(indices[0]);
      for (size_t i = 1; i + 1 < indices.size(); ++i) {
        tri_indices.push_back(v0);
        tri_indices.push_back(static_cast<std::uint32_t>(indices[i]));
        tri_indices.push_back(static_cast<std::uint32_t>(indices[i + 1]));
      }
    }
  }

  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  return true;
}

} // namespace s21