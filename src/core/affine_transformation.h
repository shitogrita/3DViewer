#pragma once

#include <cstdint>
#include <vector>

#include "matrix.h"

namespace s21{
	struct Vec3 {
		double x{};
		double y{};
		double z{};
	};

	enum class Space {
		kLocal,
		kWorldw
	};

	class AffineTransformation {
	public:
		const S21Matrix& Matrix() const noexcept;
	public:

		static S21Matrix ExpandMatrix(const S21Matrix &matrix);
		static S21Matrix ShrinkMatrix(const S21Matrix &shrinking_matrix);
		static S21Matrix Stretch(const S21Matrix& matrix, double rate);
		static S21Matrix Identity4();
		static S21Matrix Translation4(double dx, double dy, double dz);

		static S21Matrix GetRotationYMatrix(double fi);
		static S21Matrix GetRotationXMatrix(double fi);
		static S21Matrix GetRotationZMatrix(double fi);
		static S21Matrix ApplyTransformation(const S21Matrix& x, const S21Matrix& A, const Vec3& t);


	};
};
