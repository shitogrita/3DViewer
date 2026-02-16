#pragma once
#include "matrix.h"

namespace s21 {

	class Projection final {
	public:
		static S21Matrix Perspective(double vertical_fov_deg, // центральная проекция
									double aspect_ratio,
									double near_plane,
									double far_plane);

		static S21Matrix Ortho(double left, double right, // параллельная(ортографическая) проекция mvp
							   double bottom, double top,
							   double z_near, double z_far);

		static S21Matrix OrthoSymmetric(double half_width,
										double aspect,
										double z_near, double z_far);
	};

}
