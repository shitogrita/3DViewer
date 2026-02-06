#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

class S21Matrix {
private:
	int rows_, cols_ = 0;
	double **matrix_ = nullptr;

	static int NextPow2(int x);

	S21Matrix PadToSquare(int n) const;

	S21Matrix Slice(int r0, int c0, int sz) const;

	void Paste(const S21Matrix &block, int r0, int c0);

	static S21Matrix StrassenRec(const S21Matrix &A, const S21Matrix &B);

public:
	S21Matrix() {
		rows_ = 0;
		cols_ = 0;
	}

	S21Matrix(int rows, int cols);

	S21Matrix(const S21Matrix &other);

	S21Matrix(std::initializer_list<std::initializer_list<double> > list) : S21Matrix(
		list.size(), list.begin()[0].size()) {

		for (int i = 0; i < list.size(); ++i) {
			for (int j = 0; j < list.begin()[i].size(); ++j) {
				matrix_[i][j] = list.begin()[i].begin()[j];
			}
		}
	}

	S21Matrix(S21Matrix &&other);

	S21Matrix &operator=(const S21Matrix &other);

	S21Matrix &operator=(S21Matrix &&other);

	~S21Matrix();

	bool EqMatrix(const S21Matrix &other) const;

	void SumMatrix(const S21Matrix &other);

	void SubMatrix(const S21Matrix &other);

	void MulNumber(const double number);

	void MulMatrix(const S21Matrix &other);

	S21Matrix Transpose() const;

	S21Matrix CalcComplements() const;

	double Determinant() const;

	S21Matrix InverseMatrix() const;

	S21Matrix rec_minor(int temp_row, int temp_col) const;

	int get_rows() const;

	int get_cols() const;

	void set_rows(int new_rows);

	void set_cols(int new_cols);

	S21Matrix operator+(const S21Matrix &other);

	S21Matrix operator-(const S21Matrix &other);

	S21Matrix operator*(double number) const;

	S21Matrix operator*(const S21Matrix &other) const;

	bool operator==(const S21Matrix &other) const;

	S21Matrix &operator+=(const S21Matrix &other);

	S21Matrix &operator-=(const S21Matrix &other);

	S21Matrix &operator*=(const double number);

	S21Matrix &operator*=(const S21Matrix &other);

	double &operator()(int i, int j);

	double operator()(int i, int j) const;

	class Proxy {
	public:
		Proxy(double *m) : m_(m) {
		}

		double &operator[](int index) {
			return m_[index];
		}

		double &operator[](int index) const {
			return m_[index];
		}

	private:
		double *m_;
	};

	Proxy operator[](int j) { return Proxy(matrix_[j]); }

	Proxy operator[](int j) const { return Proxy(matrix_[j]); }

	void StrassenAlgorithm(const S21Matrix &other);

	double **GetMatrixArray() {
		return matrix_;
	}
};
