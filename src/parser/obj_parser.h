#pragma once

#include <string>
#include <vector>
#include <utility>
#include <charconv>
#include <fstream>
#include <algorithm>
#include <cctype>

namespace s21 {

	// Структура вершины
	struct Vertex {
		double x, y, z;
	};

	// Ребро задаётся индексами двух вершин (0‑base)
	using Edge = std::pair<unsigned int, unsigned int>;

	class ObjParser {
	public:


		// Загружает модель из .obj файла.
		// Параметры:
		//   filename  - путь к файлу
		//   vertices  - вектор, в который будут помещены вершины
		//   edges     - вектор, в который будут помещены рёбра
		// Возвращает true при успешной загрузке, false в случае ошибки.
		static bool Parse(const std::string& filename,
						  std::vector<Vertex>& vertices,
						  std::vector<Edge>& edges,
						  std::vector<std::uint32_t>& tri_indices);
	};
}

