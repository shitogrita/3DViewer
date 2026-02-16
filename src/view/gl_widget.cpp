#include "view/gl_widget.h"

#include <QKeyEvent>


namespace s21 {
  namespace {
    inline void NormalizeAngle(double& a) {
      while (a >= 360.0) a -= 360.0;
      while (a < 0.0) a += 360.0;
    }
  }

  GlWidget::GlWidget(QWidget* parent) : QOpenGLWidget(parent) {
    setFocusPolicy(Qt::StrongFocus);

    // Параметры по умолчанию: каркас + серый фон, чтобы чёрные линии были видны
    params_.fill_enabled = false;
    params_.draw_edges = true;

    params_.edge_rgb[0] = 0.f;
    params_.edge_rgb[1] = 0.f;
    params_.edge_rgb[2] = 0.f;
    params_.edge_width = 0.f;

    params_.fill_rgba[0] = 0.f; // красный
    params_.fill_rgba[1] = 0.f; // зеленый
    params_.fill_rgba[2] = 1.f; // синий!
    params_.fill_rgba[3] = 0.8f; // прозрачность

    params_.background_rgb[0] = 0.92f;
    params_.background_rgb[1] = 0.92f;
    params_.background_rgb[2] = 0.92f;

    //тестовый меш
    mesh_ = s21::test_meshes::MakeTorus(1.0f, 0.5f, 64, 32);

    connect(&timer_, &QTimer::timeout, this, [this]() {
      TickInput_();
      update();
    });
    timer_.start(16);
  }

  GlWidget::~GlWidget() {
    makeCurrent();
    render_.Destroy();
    doneCurrent();
  }

  void GlWidget::initializeGL() {
    initializeOpenGLFunctions();

    if (!render_.Initialize(this)) {
      qWarning() << "GlRender Initialize failed";
      return;
    }
    render_.UploadMesh(mesh_);
  }

  void GlWidget::resizeGL(int w, int h) {
    render_.Resize(w, h);
  }

  void GlWidget::paintGL() {
    UpdateMvp_();
    render_.Render(mvp_col_major_, params_);
  }

  void GlWidget::keyPressEvent(QKeyEvent* e) {
    if (!e->isAutoRepeat()) keys_.insert(e->key());

    if (!e->isAutoRepeat()) {
      if (e->key() == Qt::Key_P) ToggleProjection_();  // Perspective <-> Ortho
      if (e->key() == Qt::Key_Z) ToggleFill_();        // fill on/off
    }
    QOpenGLWidget::keyPressEvent(e);
  }

  void GlWidget::keyReleaseEvent(QKeyEvent* e) {
    if (!e->isAutoRepeat()) keys_.remove(e->key());
    QOpenGLWidget::keyReleaseEvent(e);
  }

  void GlWidget::TickInput_() {
    bool changed = false;
    const double rot_step = 1.0;
    const double move_step = 0.05;

    if (keys_.contains(Qt::Key_I)) { ax_ += rot_step; changed = true; }
    if (keys_.contains(Qt::Key_K)) { ax_ -= rot_step; changed = true; }
    if (keys_.contains(Qt::Key_J)) { ay_ += rot_step; changed = true; }
    if (keys_.contains(Qt::Key_L)) { ay_ -= rot_step; changed = true; }
    if (keys_.contains(Qt::Key_U)) { az_ += rot_step; changed = true; }
    if (keys_.contains(Qt::Key_O)) { az_ -= rot_step; changed = true; }
    const double zoom_k = 0.98;

    if (proj_mode_ == ProjectionMode::kPerspective) {
      if (keys_.contains(Qt::Key_W)) {cam_z_ -= move_step; changed = true;}
      if (keys_.contains(Qt::Key_S)) {cam_z_ += move_step; changed = true;}
    } else {
      if (keys_.contains(Qt::Key_W)) { ortho_half_w_ *= zoom_k; changed = true;}
      if (keys_.contains(Qt::Key_S)) { ortho_half_w_ /= zoom_k; changed = true;}
    }

    if (keys_.contains(Qt::Key_A)) { cam_x_ -= move_step; changed = true;}
    if (keys_.contains(Qt::Key_D)) { cam_x_ += move_step; changed = true;}
    if (keys_.contains(Qt::Key_Q)) { cam_y_ -= move_step; changed = true;}
    if (keys_.contains(Qt::Key_E)) { cam_y_ += move_step; changed = true;}



    if (cam_z_ < 0.2) cam_z_ = 0.2;

    if (changed) NormalizeAngles_();
  }

  void GlWidget::NormalizeAngles_() {
    NormalizeAngle(ax_);
    NormalizeAngle(ay_);
    NormalizeAngle(az_);
  }

  void GlWidget::ToggleProjection_() {  // тип проекции
    proj_mode_ = (proj_mode_ == ProjectionMode::kPerspective)
                     ? ProjectionMode::kOrtho
                     : ProjectionMode::kPerspective;
  }

  void GlWidget::ToggleFill_() { // заливка граней
    params_.fill_enabled = !params_.fill_enabled;
  }

  void GlWidget::UpdateMvp_() {
    const double aspect =
        (height() == 0) ? 1.0 : (static_cast<double>(width()) / static_cast<double>(height()));
    const S21Matrix Rx = s21::AffineTransformation::GetRotationXMatrix(ax_);
    const S21Matrix Ry = s21::AffineTransformation::GetRotationYMatrix(ay_);
    const S21Matrix Rz = s21::AffineTransformation::GetRotationZMatrix(az_);
    const S21Matrix M = Rz * (Ry * Rx);

    const S21Matrix V =
        s21::AffineTransformation::Translation4(-cam_x_, -cam_y_, -cam_z_);

    S21Matrix P(4, 4);
    if (proj_mode_ == ProjectionMode::kPerspective) {
      P = s21::Projection::Perspective(fov_y_deg_, aspect, 0.1, 100.0);

    } else {

      const double half_w = ortho_half_w_;
      const double half_h = ortho_half_w_ / aspect;
      P = s21::Projection::OrthoSymmetric(half_w, half_h, 0.1, 100.0);
    }

    const S21Matrix MVP = P * (V * M);

    mvp_col_major_ = s21::AffineTransformation::GetColMajor(MVP);
  }
}
