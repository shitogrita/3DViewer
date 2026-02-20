#pragma once

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLWidget>
#include <QSet>
#include <QTimer>

#include <array>

#include "core/affine_transformation.h"
#include "core/matrix.h"
#include "core/projection.h"
#include "../render/opengl/gl_render.h"


namespace s21 {

	class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
		Q_OBJECT

	   public:
		explicit GlWidget(QWidget* parent = nullptr);
		~GlWidget() override;

		// interface
		void SetBackgroundColor(const QColor& c);

		void SetProjectionType(int type);

		void SetEdgesDashed(bool on);
		void SetEdgeWidth(float w);
		void SetEdgeColor(const QColor& c);

		void SetVertexMode(int mode);
		void SetVertexSize(float s);
		void SetVertexColor(const QColor& c);

		void SetDashPeriod(float p);
		void SetDashFill(float f);

		void SetFillEnabled(bool on);
		void SetFillAlpha(float a);
		void SetFillColor(const QColor& c, float a);
		void SetRotationDegrees(double ax_deg, double ay_deg, double az_deg);
		void ResetRotation();

		struct RotationDeg {
			double x;
			double y;
			double z;
		};
		RotationDeg GetRotationDegrees() const;

	signals:
	  void RotationChanged(double ax_deg, double ay_deg, double az_deg);

	protected:
		void FitViewToMesh_();

		void initializeGL() override;
		void resizeGL(int w, int h) override;
		void paintGL() override;

		void keyPressEvent(QKeyEvent* e) override;
		void keyReleaseEvent(QKeyEvent* e) override;
	public:
		void SetFillOpaque(bool opaque);
	private:
		float fill_alpha_saved_ = 0.8f;

		signals:
  void FillEnabledChanged(bool on);


	private:
		bool render_ready_ = false;

	public:

		bool LoadModelFromObjFile(const QString& path);

		bool GetFillEnabled() const { return params_.fill_enabled; }
		bool GetFillOpaque() const { return params_.fill_rgba[3] >= 0.999f; }



		enum class ProjectionMode { kPerspective, kOrtho };

		void TickInput_();
		void NormalizeAngles_();

		void UpdateMvp_();
		void ToggleProjection_();
		void ToggleFill_();

		GlRender render_;
		GlRender::MeshData mesh_;
		GlRender::DrawParams params_;

		ProjectionMode proj_mode_ = ProjectionMode::kPerspective;

		// camera
		double cam_x_ = 0.0;
		double cam_y_ = 0.0;
		double cam_z_ = 3.0;

		double fov_y_deg_ = 60.0;

		// symmetric
		double ortho_half_w_ = 1.2;

		// degrees
		double ax_ = 20.0;
		double ay_ = 30.0;
		double az_ = 0.0;

		std::array<float, 16> mvp_col_major_{};

		QTimer timer_;
		QSet<int> keys_;

		signals:
		  void ScaleChanged(double s);

	public:
		double GetScale() const { return scale_; }

		void SetScale(double s);
		void ScaleBy(double k);
		void ResetScale();

	private:
		double scale_ = 1.0;

		signals:
		void ModelInfoChanged(const QString& file_name, int vertex_count, int edge_count);

	public:
		QString GetModelFileName() const { return model_file_name_; }
		int GetVertexCount() const { return model_vertex_count_; }
		int GetEdgeCount() const { return model_edge_count_; }


	private:
		QString model_file_name_ = "—";
		int model_vertex_count_ = 0;
		int model_edge_count_ = 0;
	};

}