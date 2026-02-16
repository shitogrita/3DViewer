#include "view/test_meshes.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <vector>

namespace s21::test_meshes {

namespace {

// ключ ребра (min,max) в uint64
inline std::uint64_t EdgeKey(std::uint32_t a, std::uint32_t b) {
  const std::uint32_t lo = std::min(a, b);
  const std::uint32_t hi = std::max(a, b);
  return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
}

inline void AddEdge(std::unordered_set<std::uint64_t>& s,
                    std::uint32_t a, std::uint32_t b) {
  s.insert(EdgeKey(a, b));
}

}  // namespace

GlRender::MeshData MakeTorus(float majorR, float minorR, int majorSeg, int minorSeg) {
  GlRender::MeshData m;

  majorSeg = std::max(3, majorSeg);
  minorSeg = std::max(3, minorSeg);

  const int vcount = majorSeg * minorSeg;
  m.vertices_xyz.reserve(static_cast<size_t>(vcount) * 3);

  const double twoPi = 6.28318530717958647692;

  // Вершины тора:
  // u — вокруг большого радиуса, v — вокруг малого.
  for (int i = 0; i < majorSeg; ++i) {
    const double u = twoPi * (static_cast<double>(i) / majorSeg);
    const double cu = std::cos(u);
    const double su = std::sin(u);

    for (int j = 0; j < minorSeg; ++j) {
      const double v = twoPi * (static_cast<double>(j) / minorSeg);
      const double cv = std::cos(v);
      const double sv = std::sin(v);

      const double x = (majorR + minorR * cv) * cu;
      const double y = (majorR + minorR * cv) * su;
      const double z = minorR * sv;

      m.vertices_xyz.push_back(static_cast<float>(x));
      m.vertices_xyz.push_back(static_cast<float>(y));
      m.vertices_xyz.push_back(static_cast<float>(z));
    }
  }

  // Треугольники (квад на сетке -> 2 треугольника)
  m.tri_indices.reserve(static_cast<size_t>(majorSeg) * minorSeg * 6);

  auto idx = [minorSeg](int i, int j) -> std::uint32_t {
    return static_cast<std::uint32_t>(i * minorSeg + j);
  };

  for (int i = 0; i < majorSeg; ++i) {
    const int in = (i + 1) % majorSeg;
    for (int j = 0; j < minorSeg; ++j) {
      const int jn = (j + 1) % minorSeg;

      const std::uint32_t a = idx(i,  j);
      const std::uint32_t b = idx(in, j);
      const std::uint32_t c = idx(in, jn);
      const std::uint32_t d = idx(i,  jn);

      // a-b-c и a-c-d
      m.tri_indices.push_back(a);
      m.tri_indices.push_back(b);
      m.tri_indices.push_back(c);

      m.tri_indices.push_back(a);
      m.tri_indices.push_back(c);
      m.tri_indices.push_back(d);
    }
  }

  // Рёбра из треугольников (уникальные)
  std::unordered_set<std::uint64_t> edges;
  edges.reserve(m.tri_indices.size());

  for (size_t k = 0; k + 2 < m.tri_indices.size(); k += 3) {
    const std::uint32_t a = m.tri_indices[k + 0];
    const std::uint32_t b = m.tri_indices[k + 1];
    const std::uint32_t c = m.tri_indices[k + 2];
    AddEdge(edges, a, b);
    AddEdge(edges, b, c);
    AddEdge(edges, c, a);
  }

  m.edge_indices.reserve(edges.size() * 2);
  for (std::uint64_t key : edges) {
    const std::uint32_t lo = static_cast<std::uint32_t>(key >> 32);
    const std::uint32_t hi = static_cast<std::uint32_t>(key & 0xFFFFFFFFu);
    m.edge_indices.push_back(lo);
    m.edge_indices.push_back(hi);
  }

  return m;
}

  // cube in cube
inline void AppendCubeEdges(std::vector<std::uint32_t>& out, std::uint32_t base) {
  // 12 рёбер куба (24 индекса)
  const std::uint32_t e[] = {
      0,1, 1,2, 2,3, 3,0,  // нижняя грань
      4,5, 5,6, 6,7, 7,4,  // верхняя грань
      0,4, 1,5, 2,6, 3,7   // стойки
  };
  out.reserve(out.size() + 24);
  for (int i = 0; i < 24; ++i) out.push_back(base + e[i]);
}

inline void AppendCubeTris(std::vector<std::uint32_t>& out, std::uint32_t base) {
  // 12 треугольников = 36 индексов (ориентация не критична для каркаса)
  const std::uint32_t t[] = {
      0,1,2, 0,2,3,
      4,6,5, 4,7,6,
      0,3,7, 0,7,4,
      1,5,6, 1,6,2,
      0,4,5, 0,5,1,
      3,2,6, 3,6,7
  };
  out.reserve(out.size() + 36);
  for (int i = 0; i < 36; ++i) out.push_back(base + t[i]);
}

 // namespace

GlRender::MeshData MakeInnerFrameCube(float outer_half, float inner_half) {
  GlRender::MeshData m;

  // Вершины куба: порядок соответствует вашему кубу (0..7)
  auto append_cube_vertices = [&](float h) {
    const float v[] = {
        -h, -h, -h,  // 0
         h, -h, -h,  // 1
         h,  h, -h,  // 2
        -h,  h, -h,  // 3
        -h, -h,  h,  // 4
         h, -h,  h,  // 5
         h,  h,  h,  // 6
        -h,  h,  h   // 7
    };
    m.vertices_xyz.insert(m.vertices_xyz.end(), v, v + 24);
  };

  m.vertices_xyz.reserve(24 * 2);
  append_cube_vertices(outer_half);  // 0..7
  append_cube_vertices(inner_half);  // 8..15

  // Рёбра: внешний куб + внутренний куб + 8 диагоналей (соответствующие вершины)
  std::vector<std::uint32_t> edges;
  edges.reserve(24 + 24 + 16);

  AppendCubeEdges(edges, /*base=*/0);   // внешний
  AppendCubeEdges(edges, /*base=*/8);   // внутренний

  // Диагональные распорки: i внешняя -> i внутренняя
  for (std::uint32_t i = 0; i < 8; ++i) {
    edges.push_back(i);
    edges.push_back(8 + i);
  }

  m.edge_indices = std::move(edges);

  // Треугольники (если захотите включать fill). Можно:
  // 1) заливать только внешний и внутренний куб,
  // 2) либо ещё добавлять “стенки” между ними (это уже другой объект).
  // Для теста достаточно (1).
  m.tri_indices.reserve(36 + 36);
  AppendCubeTris(m.tri_indices, 0);
  AppendCubeTris(m.tri_indices, 8);

  return m;

}

}
