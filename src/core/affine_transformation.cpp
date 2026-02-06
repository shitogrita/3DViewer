#include "affine_transformation.h"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace s21 {
	S21Matrix AffineTransformation::Identity4() {
		S21Matrix m(4,4);
		for (int i = 0; i < 4 ; ++i) {
			for (int j = 0; j < 4; ++ j) {
				m (i,j) = (i == j) ? 1.0 : 0.0;
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
			for (int j = 0 ; j < expanding_matrix.get_cols(); j ++) {
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
			for (int j = 0 ; j < shrinking_matrix.get_cols(); j ++) {
				matrix[i][j] = shrinking_matrix[i][j];
			}
		}
		return matrix;
	}


	S21Matrix AffineTransformation::ApplyTransformation(const S21Matrix& x, const S21Matrix& A, const Vec3& t) {

		S21Matrix T(4, 4);

		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				T[i][j] = A[i][j];
			}
		}

		T[0][3] = t.x;
		T[1][3] = t.y;
		T[2][3] = t.z;
		T[3][0] = 0.0;
		T[3][1] = 0.0;
		T[3][2] = 0.0;
		T[3][3] = 1.0;

		return ShrinkMatrix(T * ExpandMatrix(x));
	}

	S21Matrix AffineTransformation::GetRotationYMatrix(double fi) {
		fi = fi *  M_PI / 180;
		return {
			{std::cos(fi), 0.0, std::sin(fi), 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{-std::sin(fi), 0.0, std::cos(fi), 0.0},
			{0.0, 0.0, 0.0, 1.0}
		};
	}

	S21Matrix AffineTransformation::GetRotationXMatrix(double fi) {
		fi = fi *  M_PI / 180;
		return {
				{1.0, 0.0, 0.0 , 0.0},
				{0.0, std::cos(fi), -std::sin(fi), 0.0},
				{0.0, std::sin(fi),  std::cos(fi), 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
	}

	S21Matrix AffineTransformation::GetRotationZMatrix(double fi) {
		fi = fi *  M_PI / 180;
		return {
				{std::cos(fi), -std::sin(fi), 0.0, 0.0},
				{std::sin(fi), std::cos(fi), 0.0, 0.0},
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
		};
	}


	S21Matrix AffineTransformation::Stretch(const S21Matrix& matrix, double rate) {
		auto result = matrix * rate;
		return result;
	}
}

