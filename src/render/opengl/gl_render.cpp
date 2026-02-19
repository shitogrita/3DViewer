#include "gl_render.h"

#include <iostream>

using s21::GlRender;

GlRender::~GlRender() { Destroy(); }

namespace s21 {

	bool GlRender::Initialize(QOpenGLFunctions_3_3_Core* f,
							const QString& vertex_shader,
							const QString& fragment_shader) {

		if (QOpenGLContext::currentContext() == nullptr) {
			exit(3);
		}
		f_ = f;
		if (!f_) return false;

		f_->glEnable(GL_DEPTH_TEST);
		f_->glDepthFunc(GL_LESS);

		f_->glEnable(GL_PROGRAM_POINT_SIZE);
		f_->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		QString log;
		if (!GlShaderLoader::load(program_, vertex_shader, fragment_shader, &log)) {
			qWarning().noquote() << "Shader loader failed.\n"
								 << "VS path =" << vertex_shader << "\n"
								 << "FS path =" << fragment_shader << "\n"
								 << log;
			return false;
		}

		program_.bind();

		loc_mvp_ = program_.uniformLocation("uMVP"); //точное определение в памяти
		loc_color_ = program_.uniformLocation("uColor");

		loc_use_dash_    = program_.uniformLocation("uUseDash");
		loc_dash_period_ = program_.uniformLocation("uDashPeriod");
		loc_dash_fill_   = program_.uniformLocation("uDashFill");

		loc_point_mode_  = program_.uniformLocation("uPointMode");
		loc_point_soft_  = program_.uniformLocation("uPointSoft");
		loc_point_size_  = program_.uniformLocation("uPointSize");
		program_.release();

		//GPU
		f_->glGenVertexArrays(1, &vao_);
		f_->glGenBuffers(1, &vbo_);
		f_->glGenBuffers(1, &ebo_tri_);
		f_->glGenBuffers(1, &ebo_edge_);
		f_->glGenVertexArrays(1, &vao_edge_);
		f_->glGenBuffers(1, &vbo_edge_);
		return true;
	}

	void GlRender::UploadMesh(const MeshData& mesh) {
		if (!f_) return;

		tri_index_count_ = static_cast<int>(mesh.tri_indices.size());
		edge_index_count_ = static_cast<int>(mesh.edge_indices.size());
		vertex_count_ = static_cast<int>(mesh.vertices_xyz.size() / 3);

		f_->glBindVertexArray(vao_);

		f_->glBindBuffer(GL_ARRAY_BUFFER, vbo_);
		f_->glBufferData(GL_ARRAY_BUFFER,
						 static_cast<GLsizeiptr>(mesh.vertices_xyz.size() * sizeof(float)),
						 mesh.vertices_xyz.data(),
						 GL_STATIC_DRAW);

		f_->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
		f_->glEnableVertexAttribArray(0);

		f_->glDisableVertexAttribArray(1);

		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_tri_);
		f_->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
						static_cast<GLsizeiptr>(mesh.tri_indices.size() * sizeof(std::uint32_t)),
						mesh.tri_indices.data(),
						GL_STATIC_DRAW);

		f_->glBindVertexArray(0);

		std::vector<EdgeVertex> edge_verts;
		edge_verts.reserve(mesh.edge_indices.size());

		 auto get_pos = [&](std::uint32_t idx) { // индексы -> ребра
			 const float* p = &mesh.vertices_xyz[idx*3];
		 	return std::array<float, 3>{p[0],p[1],p[2]};
		 };
		for (size_t k = 0; k + 1 < mesh.edge_indices.size(); k += 2) {
			std::uint32_t i0 = mesh.edge_indices[k];
			std::uint32_t i1 = mesh.edge_indices[k+1];

			auto p0 = get_pos(i0);
			auto p1 = get_pos(i1);

			float dx = p1[0] - p0[0];
			float dy = p1[1] - p0[1];
			float dz = p1[2] - p0[2];
			float L  = std::sqrt(dx * dx + dy * dy + dz * dz); // длина ребра

			edge_verts.push_back({p0[0], p0[1], p0[2], 0.0f});
			edge_verts.push_back({p1[0], p1[1], p1[2], L});
		}

		edge_vertex_count_ = static_cast<int>(edge_verts.size());

		// upload GPU
		f_->glBindVertexArray(vao_edge_);
		f_->glBindBuffer(GL_ARRAY_BUFFER, vbo_edge_);
		f_->glBufferData(GL_ARRAY_BUFFER,
						 static_cast<GLsizeiptr>(edge_verts.size() * sizeof(EdgeVertex)),
						 edge_verts.data(),
						 GL_STATIC_DRAW);

		// loc = 0 vec3 aPos
		f_->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(EdgeVertex), (void*)0);
		f_->glEnableVertexAttribArray(0);

		// loc = 1: float aEdgeT
		f_->glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(EdgeVertex),
								  (void*)(3 * sizeof(float)));
		f_->glEnableVertexAttribArray(1);

		f_->glBindVertexArray(0);

		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_edge_);
		f_->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
						 static_cast<GLsizeiptr>(mesh.edge_indices.size() * sizeof(std::uint32_t)),
						 mesh.edge_indices.data(),
						 GL_STATIC_DRAW);
		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		qDebug() << "tri_index_count_ =" << tri_index_count_
		 << "edge_index_count_ =" << edge_index_count_
		 << "edge_vertex_count_ =" << edge_vertex_count_
		 << "vertex_count_ =" << vertex_count_;
	}

	void GlRender::Resize(int w, int h) {
		if (!f_) return;
		f_->glViewport(0, 0, w, h);
	}

	void GlRender::Render(const std::array<float, 16>& mvp_col_major,
					  const DrawParams& params) {
		if (!f_) return;

		// фон
		f_->glClearColor(params.background_rgb[0],
						 params.background_rgb[1],
						 params.background_rgb[2],
						 1.f);
		//f_->glClearColor(0.8f, 0.8f, 0.8f, 1.f); - серый фон
		f_->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		program_.bind();

		f_->glUniformMatrix4fv(loc_mvp_, 1, GL_FALSE, mvp_col_major.data());

		if (loc_use_dash_ != -1)   f_->glUniform1i(loc_use_dash_, 0); // нейтральные значения режимов
		if (loc_point_mode_ != -1) f_->glUniform1i(loc_point_mode_, 0);
		if (loc_point_size_ != -1) f_->glUniform1f(loc_point_size_, 1.0f);
		if (loc_point_soft_ != -1) f_->glUniform1f(loc_point_soft_, 0.0f);

		// грани
		f_->glBindVertexArray(vao_);

		if (params.fill_enabled && tri_index_count_ > 0) {
			const float a = params.fill_rgba[3];
			const bool use_blend = params.transparent || (a < 1.f);

			if (use_blend) {
				f_->glEnable(GL_BLEND);
				f_->glDepthMask(GL_FALSE);
			} else {
				f_->glDisable(GL_BLEND);
				f_->glDepthMask(GL_TRUE);
			}

			f_->glEnable(GL_POLYGON_OFFSET_FILL);
			f_->glPolygonOffset(1.f, 1.f);

			f_->glUniform4f(loc_color_,
							params.fill_rgba[0], params.fill_rgba[1],
							params.fill_rgba[2], params.fill_rgba[3]);

			f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_tri_);
			f_->glDrawElements(GL_TRIANGLES, tri_index_count_, GL_UNSIGNED_INT, nullptr);

			f_->glDisable(GL_POLYGON_OFFSET_FILL);

			f_->glDepthMask(GL_TRUE);
			f_->glDisable(GL_BLEND);
		}

		// ребра
		if (params.draw_edges && edge_vertex_count_ > 0) {
			f_->glUniform4f(loc_color_,
							params.edge_rgb[0], params.edge_rgb[1], params.edge_rgb[2], 1.f);

			f_->glLineWidth(params.edge_width);

			if (loc_use_dash_ != -1)    f_->glUniform1i(loc_use_dash_, params.edges_dashed ? 1 : 0);
			if (loc_dash_period_ != -1) f_->glUniform1f(loc_dash_period_, params.dash_period);
			if (loc_dash_fill_ != -1)   f_->glUniform1f(loc_dash_fill_, params.dash_fill);

			f_->glBindVertexArray(vao_edge_);
			f_->glDrawArrays(GL_LINES, 0, edge_vertex_count_);
			f_->glBindVertexArray(0);

			if (loc_use_dash_ != -1) f_->glUniform1i(loc_use_dash_, 0);
		}

		// вершины
		if (params.vertex_mode != 0 && vertex_count_ > 0) {
			f_->glEnable(GL_BLEND); // мягкий край

			f_->glUniform4f(loc_color_,
							params.vertex_rgb[0], params.vertex_rgb[1], params.vertex_rgb[2], 1.f);

			if (loc_point_mode_ != -1) f_->glUniform1i(loc_point_mode_, params.vertex_mode); // 1 круг, 2 квадрат
			if (loc_point_size_ != -1) f_->glUniform1f(loc_point_size_, params.vertex_size);
			if (loc_point_soft_ != -1) f_->glUniform1f(loc_point_soft_, 0.08f);

			f_->glBindVertexArray(vao_);
			f_->glDrawArrays(GL_POINTS, 0, vertex_count_);
			f_->glBindVertexArray(0);

			if (loc_point_mode_ != -1) f_->glUniform1i(loc_point_mode_, 0);
			f_->glDisable(GL_BLEND);
		}
		f_->glBindVertexArray(0);
		program_.release();
	}

	void GlRender::Destroy() {
		if (!f_) return;

		if (vbo_edge_) f_->glDeleteBuffers(1, &vbo_edge_);
		if (vao_edge_) f_->glDeleteVertexArrays(1, &vao_edge_);

		if (ebo_edge_) f_->glDeleteBuffers(1, &ebo_edge_);
		if (ebo_tri_) f_->glDeleteBuffers(1, &ebo_tri_);
		if (vbo_) f_->glDeleteBuffers(1, &vbo_);
		if (vao_) f_->glDeleteVertexArrays(1, &vao_);

		vbo_edge_ = vao_edge_ = 0;
		ebo_edge_ = ebo_tri_ = vbo_ = vao_ = 0;

		tri_index_count_ = 0;
		edge_index_count_ = 0;
		edge_vertex_count_ = 0;
		vertex_count_ = 0;
	}

}

