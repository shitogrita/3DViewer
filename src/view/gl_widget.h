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

#include "view/test_meshes.h"


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
		void initializeGL() override;
		void resizeGL(int w, int h) override;
		void paintGL() override;

		void keyPressEvent(QKeyEvent* e) override;
		void keyReleaseEvent(QKeyEvent* e) override;

	private:
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
	};

}