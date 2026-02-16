#pragma once

#include <QString>

class QOpenGLShaderProgram;

namespace s21 {

	class GlShaderLoader final {
	public:
		static bool load(QOpenGLShaderProgram& program,
						 const QString& vertex_path,
						 const QString& fragment_path,
						 QString* out_error_log);
	};

}
