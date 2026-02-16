#include "projection.h"

#include <cmath>

namespace s21 {

  S21Matrix Projection::Perspective(double fov_y_deg, double aspect,
                                  double z_near, double z_far) {
    const double fov_rad = fov_y_deg * M_PI / 180.0;
    const double f = 1.0 / std::tan(fov_rad * 0.5);

    S21Matrix p(4, 4);

    p(0, 0) = f / aspect;
    p(0, 1) = 0.0;
    p(0, 2) = 0.0;
    p(0, 3) = 0.0;

    p(1, 0) = 0.0;
    p(1, 1) = f;
    p(1, 2) = 0.0;
    p(1, 3) = 0.0;

    p(2, 0) = 0.0;
    p(2, 1) = 0.0;
    p(2, 2) = (z_far + z_near) / (z_near - z_far);
    p(2, 3) = (2.0 * z_far * z_near) / (z_near - z_far);

    p(3, 0) = 0.0;
    p(3, 1) = 0.0;
    p(3, 2) = -1.0;
    p(3, 3) = 0.0;

    return p;
  }

  S21Matrix Projection::Ortho(double left, double right,
                              double bottom, double top,
                              double z_near, double z_far) {
    S21Matrix o(4, 4);

    o(0, 0) =  2.0 / (right - left);
    o(0, 1) =  0.0;
    o(0, 2) =  0.0;
    o(0, 3) = -(right + left) / (right - left);

    o(1, 0) =  0.0;
    o(1, 1) =  2.0 / (top - bottom);
    o(1, 2) =  0.0;
    o(1, 3) = -(top + bottom) / (top - bottom);

    o(2, 0) =  0.0;
    o(2, 1) =  0.0;
    o(2, 2) = -2.0 / (z_far - z_near);
    o(2, 3) = -(z_far + z_near) / (z_far - z_near);

    o(3, 0) = 0.0;
    o(3, 1) = 0.0;
    o(3, 2) = 0.0;
    o(3, 3) = 1.0;

    return o;
  }

  S21Matrix Projection::OrthoSymmetric(double half_w, double half_h,
                                       double z_near, double z_far) {
    return Ortho(-half_w, +half_w, -half_h, +half_h, z_near, z_far);
  }

}
