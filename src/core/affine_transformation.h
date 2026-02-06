#pragma once

#include "matrix.h"

namespace s21{
	struct Vec3 {
		double x{};
		double y{};
		double z{};
	};

	class AffineTransformation {
	public:
		static S21Matrix Identity4();
		static S21Matrix Identity3();
		static S21Matrix RotateBy(const S21Matrix& matrix, double OX_degree, double OY_degree, double OZ_degree);
		static S21Matrix MoveObject(const S21Matrix& matrix, double x, double y, double z);
		static S21Matrix Stretch(const S21Matrix& matrix, double rate);

	private:
		static S21Matrix ExpandMatrix(const S21Matrix &expanding_matrix);
		static S21Matrix ShrinkMatrix(const S21Matrix &shrinking_matrix);

		static S21Matrix Translation4(double dx, double dy, double dz);

		static S21Matrix GetRotationYMatrix(double fi);
		static S21Matrix GetRotationXMatrix(double fi);
		static S21Matrix GetRotationZMatrix(double fi);
		static S21Matrix ApplyTransformation(const S21Matrix& x, const S21Matrix& A, const Vec3& t);


		static s21::Vec3 ComputeBBoxCenter(const S21Matrix& v) {
			// v: 3xN
			if (v.get_rows() != 3) {
				throw std::invalid_argument("ComputeBBoxCenter: expected 3xN matrix");
			}

			const int n = v.get_cols();

			double min_x = v[0][0], max_x = v[0][0];
			double min_y = v[1][0], max_y = v[1][0];
			double min_z = v[2][0], max_z = v[2][0];

			for (int j = 1; j < n; ++j) {
				const double x = v[0][j];
				const double y = v[1][j];
				const double z = v[2][j];

				if (x < min_x) min_x = x;
				if (x > max_x) max_x = x;

				if (y < min_y) min_y = y;
				if (y > max_y) max_y = y;

				if (z < min_z) min_z = z;
				if (z > max_z) max_z = z;
			}

			return s21::Vec3{
				(min_x + max_x) * 0.5,
				(min_y + max_y) * 0.5,
				(min_z + max_z) * 0.5
			};
		}

	public:

		static void PrintMatrix(const S21Matrix& result, bool should_transpose = false) {

			auto z = [](double v) {
				return (std::abs(v) < 1e-12) ? 0.0 : v;
			};
			auto x = should_transpose ? result.Transpose() : result;
			for (int i = 0; i < x.get_rows(); ++i) {
				for (int j = 0; j < x.get_cols(); ++j) {
					std::cout << z(x[i][j]) << '\t';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;

		}
	};
};
