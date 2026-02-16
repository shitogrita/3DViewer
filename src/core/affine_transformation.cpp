#include "affine_transformation.h"

#include <cmath>

namespace s21 {
	S21Matrix AffineTransformation::Identity4() {
		S21Matrix m(4, 4);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				m(i, j) = (i == j) ? 1.0 : 0.0;
			}
		}
		return m;
	}

	S21Matrix AffineTransformation::Translation4(double dx, double dy, double dz) {
		S21Matrix t = Identity4();
		t(0, 3) = dx;
		t(1, 3) = dy;
		t(2, 3) = dz;
		return t;
	}

	S21Matrix AffineTransformation::ExpandMatrix(const S21Matrix &expanding_matrix) {
		S21Matrix matrix(4, expanding_matrix.get_cols());
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < expanding_matrix.get_cols(); j++) {
				matrix[i][j] = expanding_matrix[i][j];
			}
		}
		for (int i = 0; i < expanding_matrix.get_cols(); i++) {
			matrix[3][i] = 1.0;
		}

		return matrix;
	}

	S21Matrix AffineTransformation::ShrinkMatrix(const S21Matrix &shrinking_matrix) {
		S21Matrix matrix(3, shrinking_matrix.get_cols());
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < shrinking_matrix.get_cols(); j++) {
				matrix[i][j] = shrinking_matrix[i][j];
			}
		}
		return matrix;
	}

	S21Matrix AffineTransformation::BuildAffine4(const S21Matrix& A, const Vec3& t) {
		S21Matrix T(4, 4);
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				T(i, j) = A(i, j);
			}
		}
		T(0, 3) = t.x;
		T(1, 3) = t.y;
		T(2, 3) = t.z;
		T(3, 0) = 0.0; T(3, 1) = 0.0; T(3, 2) = 0.0; T(3, 3) = 1.0;

		return T;
	}

	S21Matrix AffineTransformation::ApplyTransformation(const S21Matrix& x,
													const S21Matrix& A,
													const Vec3& t) {
		const S21Matrix T = BuildAffine4(A, t);
		return ShrinkMatrix(T * ExpandMatrix(x));
	}

	S21Matrix AffineTransformation::GetRotationYMatrix(double fi) {
		fi = fi * M_PI / 180.0;
		return {
			{std::cos(fi), 0.0, std::sin(fi), 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{-std::sin(fi), 0.0, std::cos(fi), 0.0},
			{0.0, 0.0, 0.0, 1.0}
		};
	}

	S21Matrix AffineTransformation::GetRotationXMatrix(double fi) {
		fi = fi * M_PI / 180.0;
		return {
			{1.0, 0.0, 0.0, 0.0},
			{0.0, std::cos(fi), -std::sin(fi), 0.0},
			{0.0, std::sin(fi), std::cos(fi), 0.0},
			{0.0, 0.0, 0.0, 1.0}
		};
	}

	S21Matrix AffineTransformation::GetRotationZMatrix(double fi) {
		fi = fi * M_PI / 180.0;
		return {
			{std::cos(fi), -std::sin(fi), 0.0, 0.0},
			{std::sin(fi), std::cos(fi), 0.0, 0.0},
			{0.0, 0.0, 1.0, 0.0},
			{0.0, 0.0, 0.0, 1.0}
		};
	}

	S21Matrix AffineTransformation::Stretch(const S21Matrix &matrix, double rate) {
		auto result = matrix * rate;
		return result;
	}

	S21Matrix AffineTransformation::RotateBy(const S21Matrix &matrix,
											double OX_degree,
											double OY_degree,
											double OZ_degree) {
		const Vec3 c = ComputeBBoxCenter(matrix);

		//auto moved0 = ApplyTransformation(matrix, Identity4(), Vec3{-c.x, -c.y, -c.z});
		//auto r1 = ApplyTransformation(moved0, GetRotationXMatrix(OX_degree), {});
		//auto r2 = ApplyTransformation(r1, GetRotationYMatrix(OY_degree), {});
		//auto r3 = ApplyTransformation(r2, GetRotationZMatrix(OZ_degree), {});
		//auto moved1 = ApplyTransformation(r3, Identity4(), Vec3{c.x, c.y, c.z});
		//return moved1;

		// M = T(c) * Rz * Ry * Rx * T(-c)
		const S21Matrix T0 = Translation4(-c.x, -c.y, -c.z);
		const S21Matrix Rx = GetRotationXMatrix(OX_degree);
		const S21Matrix Ry = GetRotationYMatrix(OY_degree);
		const S21Matrix Rz = GetRotationZMatrix(OZ_degree);
		const S21Matrix T1 = Translation4(c.x, c.y, c.z);

		const S21Matrix M = T1 * (Rz * (Ry * (Rx * T0)));

		return ShrinkMatrix(M * ExpandMatrix(matrix));

	}

	S21Matrix AffineTransformation::ModelRotateAroundCenter4(const S21Matrix& vertices3xN,
														 double ox_deg, double oy_deg, double oz_deg) {
		const Vec3 c = ComputeBBoxCenter(vertices3xN);

		const S21Matrix T0 = Translation4(-c.x, -c.y, -c.z);
		const S21Matrix Rx = GetRotationXMatrix(ox_deg);
		const S21Matrix Ry = GetRotationYMatrix(oy_deg);
		const S21Matrix Rz = GetRotationZMatrix(oz_deg);
		const S21Matrix T1 = Translation4(c.x, c.y, c.z);

		return T1 * (Rz * (Ry * (Rx * T0)));
	}

	S21Matrix AffineTransformation::Identity3() {
		S21Matrix m(3, 3);
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				m(i, j) = (i == j) ? 1.0 : 0.0;
			}
		}
		return m;
	}

	S21Matrix AffineTransformation::MoveObject(const S21Matrix& matrix, double x, double y, double z) {
		return  ApplyTransformation(matrix, Identity3() , {x, y, z} );
	}

	std::array<float, 16> AffineTransformation::GetColMajor(const S21Matrix& matrix) {
		std::array<float, 16> out{};
		for (int r = 0; r < 4; ++r) {
			for (int c = 0; c < 4; ++c) {
				out[c * 4 + r] = static_cast<float>(matrix(r, c));
			}
		}
		return out;
	}
}
