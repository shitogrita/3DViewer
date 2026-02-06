#include <core/matrix.h>
#include <core/affine_transformation.h>

int main() {
	auto z = [](double v) {
		return (std::abs(v) < 1e-12) ? 0.0 : v;
	};

	S21Matrix matrix = {
		{0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0}, // x
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, // y
		{0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0} // z
	};


	auto result = s21::AffineTransformation::ApplyTransformation(
		matrix, s21::AffineTransformation::Identity4() * 5.0, {5, 5, 5});
	for (int i = 0; i < result.get_rows(); ++i) {
		for (int j = 0; j < result.get_cols(); ++j) {
			std::cout << z(result[i][j]) << " ";
		}
		std::cout << std::endl;
	}
}
