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
		program_.release();

		//GPU
		f_->glGenVertexArrays(1, &vao_);
		f_->glGenBuffers(1, &vbo_);
		f_->glGenBuffers(1, &ebo_tri_);
		f_->glGenBuffers(1, &ebo_edge_);
		return true;
	}

	void GlRender::UploadMesh(const MeshData& mesh) {
		if (!f_) return;

		tri_index_count_ = static_cast<int>(mesh.tri_indices.size());
		edge_index_count_ = static_cast<int>(mesh.edge_indices.size());

		f_->glBindVertexArray(vao_);

		f_->glBindBuffer(GL_ARRAY_BUFFER, vbo_);
		f_->glBufferData(GL_ARRAY_BUFFER,
						 static_cast<GLsizeiptr>(mesh.vertices_xyz.size() * sizeof(float)),
						 mesh.vertices_xyz.data(),
						 GL_STATIC_DRAW);

		f_->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(float), (void*)0);
		f_->glEnableVertexAttribArray(0);

		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_tri_);
		f_->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
						static_cast<GLsizeiptr>(mesh.tri_indices.size() * sizeof(std::uint32_t)),
						mesh.tri_indices.data(),
						GL_STATIC_DRAW);

		f_->glBindVertexArray(0);

		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_edge_);
		f_->glBufferData(GL_ELEMENT_ARRAY_BUFFER,
						 static_cast<GLsizeiptr>(mesh.edge_indices.size() * sizeof(std::uint32_t)),
						 mesh.edge_indices.data(),
						 GL_STATIC_DRAW);
		f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
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
		f_->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		program_.bind();
		f_->glUniformMatrix4fv(loc_mvp_, 1, GL_FALSE, mvp_col_major.data());

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

		if (params.draw_edges && edge_index_count_ > 0) {
			f_->glUniform4f(loc_color_,
							params.edge_rgb[0], params.edge_rgb[1], params.edge_rgb[2], 1.f);

			f_->glLineWidth(params.edge_width);

			f_->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_edge_);
			f_->glDrawElements(GL_LINES, edge_index_count_, GL_UNSIGNED_INT, nullptr);
		}

		f_->glBindVertexArray(0);
		program_.release();
	}


	void GlRender::Destroy() {
		if (!f_) return;

		if (ebo_edge_) f_->glDeleteBuffers(1, &ebo_edge_);
		if (ebo_tri_)  f_->glDeleteBuffers(1, &ebo_tri_);
		if (vbo_)      f_->glDeleteBuffers(1, &vbo_);
		if (vao_)      f_->glDeleteVertexArrays(1, &vao_);

		ebo_edge_ = ebo_tri_ = vbo_ = vao_ = 0;
		tri_index_count_ = 0;
		edge_index_count_ = 0;
	}

}

