#pragma once

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QString>

#include <array>
#include <cstdint>
#include <vector>

#include "gl_shader_loader.h"

namespace s21 {
	class GlRender {


	public:
		struct MeshData {
			std::vector<float> vertices_xyz;
			std::vector<std::uint32_t> tri_indices; // заливка граней
			std::vector<std::uint32_t> edge_indices; // индексы ребер
		};

		static s21::GlRender::MeshData MakeWireCube() {
			s21::GlRender::MeshData m;

			// 8 вершин куба
			const float v[] = {
				-0.5f,-0.5f,-0.5f,  // 0
				 0.5f,-0.5f,-0.5f,  // 1
				 0.5f, 0.5f,-0.5f,  // 2
				-0.5f, 0.5f,-0.5f,  // 3
				-0.5f,-0.5f, 0.5f,  // 4
				 0.5f,-0.5f, 0.5f,  // 5
				 0.5f, 0.5f, 0.5f,  // 6
				-0.5f, 0.5f, 0.5f   // 7
			  };
			m.vertices_xyz.assign(v, v + 8 * 3);

			// 12 рёбер -> 24 индекса (GL_LINES)
			const std::uint32_t e[] = {
				0,1, 1,2, 2,3, 3,0,
				4,5, 5,6, 6,7, 7,4,
				0,4, 1,5, 2,6, 3,7
			  };
			m.edge_indices.assign(e, e + 24);

			return m;
		}

		struct DrawParams {
			bool fill_enabled = false; // заливка граней
			bool draw_edges = true; // ребра

			float fill_rgba[4] = {0.f, 0.f, 1.f, 1.f}; //цвет заливки граней
			bool transparent = false;

			float edge_rgb[3] = {0.f, 0.f, 0.f}; //цвет ребер
			float edge_width = 2.f; // толщина заливки

			// bool draw_vertices = false;   // показывать вершины
			// int vertex_style = 0;         // 0=none,1=circle,2=square (реализуется отдельным draw pass)
			// float vertex_size = 5.f;      // размер вершин
			// float vertex_rgb[3] = {...};  // цвет вершин
			float background_rgb[3] = {1.f, 1.f, 1.f};  // цвет фона
		};

		GlRender() = default;
		~GlRender();

		bool Initialize(QOpenGLFunctions_3_3_Core* f,
					const QString& vertex_shader = "resources/shaders/basic.vert",
					const QString& fragment_shader = "resources/shaders/basic.frag");


		// парсера:
		//  - парсер .obj (ветка Front/Model) создаёт ModelMesh,
		//  - Controller конвертирует ModelMesh -> MeshData,
		//  - затем вызывает gl_render.UploadMesh(meshData).
		void UploadMesh(const MeshData& mesh);

		void Resize(int w, int h);

		void Render(const std::array<float, 16>& mvp_col_major, const DrawParams& params);

		void Destroy();
	private:
		QOpenGLFunctions_3_3_Core* f_ = nullptr;

		QOpenGLShaderProgram program_; // вообще после последнего апдейта можно удалить
		int loc_mvp_ = -1; // location for uMVP
		int loc_color_ = -1; // location for color (for 1)

		unsigned vao_ = 0; // формат вершины
		unsigned vbo_ = 0; // массив вершин
		unsigned ebo_tri_ = 0; // заливка треугольников
		unsigned ebo_edge_ = 0; // каркас

		int tri_index_count_ = 0;
		int edge_index_count_ = 0;
	};
}
