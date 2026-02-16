#include "gl_shader_loader.h"

#include <QFile>
#include <QOpenGLShader>
#include <QOpenGLShaderProgram>
#include <QTextStream>

namespace s21 {
	namespace {
		struct FileReadResult {
			bool ok = false;
			QString text;
			QString reason;
		};
		FileReadResult ReadTextFile(const QString& path) {
			FileReadResult r;

			QFile f(path);
			if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
				r.ok = false;
				r.reason = f.errorString();
				return r;
			}

			QTextStream in(&f);
			r.text = in.readAll();
			r.ok = true;
			return r;
		}
		QString MakeFileError(const QString& kind,
							 const QString& path,
							 const QString& reason) {

			return QString("Shader load failed: cannot read shader file. kind=%1, path=%2, reason=%3")
				.arg(kind, path, reason);
		}
		QString MakeCompileError(const QString& kind, const QString& log) {
			return QString("Shader load failed: %1 shader compile error:\n%2").arg(kind, log);
		}
		QString MakeLinkError(const QString& log) {
			return QString("Shader load failed: program link error:\n%1").arg(log);
		}
	}

	bool GlShaderLoader::load(QOpenGLShaderProgram& program,
								const QString& vertex_path,
								const QString& fragment_path,
								QString* out_error_log) {
		program.removeAllShaders();

		const auto vs_res = ReadTextFile(vertex_path);
		if (!vs_res.ok) {
			if (out_error_log) *out_error_log = MakeFileError("vertex", vertex_path, vs_res.reason);
			return false;
		}
		const auto fs_res = ReadTextFile(fragment_path);
		if (!fs_res.ok) {
			if (out_error_log) *out_error_log = MakeFileError("fragment", fragment_path, vs_res.reason);
			return false;
		}
		if (!program.addShaderFromSourceCode(QOpenGLShader::Vertex, vs_res.text)) {
			if (out_error_log) *out_error_log = MakeCompileError("vertex", program.log());
			return false;
		}

		if (!program.addShaderFromSourceCode(QOpenGLShader::Fragment, fs_res.text)) {
			if (out_error_log) *out_error_log = MakeCompileError("fragment", program.log());
			return false;
		}

		if (!program.link()) {
			if (out_error_log) *out_error_log = MakeLinkError(program.log());
			return false;
		}

		return true;
	}
}