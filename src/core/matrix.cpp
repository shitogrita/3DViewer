#include "matrix.h"

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("Incorrect matrix");
  }
  if (rows == 0 || cols == 0) {
    matrix_ = nullptr;
  } else {
    matrix_ = new double*[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_]();
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
  }
  for (int i = 0; i < rows_; ++i) {
    std::copy(other.matrix_[i], other.matrix_[i] + cols_, matrix_[i]);
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  // other.~S21Matrix();
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) return *this;

  // очистка
  if (matrix_) {
    for (int i = 0; i < rows_; ++i) delete[] matrix_[i];
    delete[] matrix_;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  if (rows_ == 0 || cols_ == 0) {
    matrix_ = nullptr;
    return *this;
  }

  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_];
    std::copy(other.matrix_[i], other.matrix_[i] + cols_, matrix_[i]);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    if (matrix_) {
      for (int i = 0; i < rows_; ++i) {
        delete[] matrix_[i];
      }
      delete[] matrix_;
    }
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = std::move(other.matrix_);
    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i < rows_; ++i) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(other.matrix_[i][j] - matrix_[i][j]) > 1e-7) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double number) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= number;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.get_rows()) {
    throw std::logic_error(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }
  S21Matrix result{rows_, other.get_cols()};
  for (int i = 0; i < result.get_rows(); ++i) {
    for (int j = 0; j < result.get_cols(); ++j) {
      for (int k = 0; k < cols_; ++k) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result = S21Matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

S21Matrix S21Matrix::rec_minor(int temp_row, int temp_col) const {
  if (temp_row < 0 || temp_row >= rows_ || temp_col < 0 || temp_col >= cols_) {
    throw std::out_of_range("Excluded row or column out of range.");
  }
  S21Matrix minor(rows_ - 1, cols_ - 1);
  int minor_i = 0;
  for (int i = 0; i < rows_; ++i) {
    if (i == temp_row) {
      continue;
    }
      int minor_j = 0;
      for (int j = 0; j < cols_; ++j) {
        if (j == temp_col) {
          continue;
        }
        minor.matrix_[minor_i][minor_j] = matrix_[i][j];
        minor_j++;
      }
      minor_i++;
    }
  return minor;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("The matrix is not square.");
  }
  S21Matrix result{rows_, cols_};
  for (int i = 0; i < result.get_rows(); ++i) {
    for (int j = 0; j < result.get_cols(); ++j) {
      S21Matrix minor = rec_minor(i, j);
      double det = minor.Determinant();
      result.matrix_[i][j] = det * ((i + j) % 2 == 0 ? 1 : -1);
    }
  }
  return result;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("The matrix is not square.");
  }
  if (rows_ == 1) {
    return matrix_[0][0];
  }
  if (rows_ == 2) {
    return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  }
  double det = 0.0;
  for (int g = 0; g < rows_; g++) {
    S21Matrix minor{rows_ - 1, rows_ - 1};
    for (int i = 1; i < rows_; i++) {
      int sub_col = 0;
      for (int j = 0; j < rows_; j++) {
        if (j == g) {
          continue;
        }
        minor.matrix_[i - 1][sub_col] = matrix_[i][j];
        ++sub_col;
      }
    }
    det += (g % 2 == 0 ? 1 : -1) * matrix_[0][g] * minor.Determinant();
  }
  return det;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::logic_error("Matrix determinant is 0.");
  }
  double det = Determinant();
  if (std::abs(det) < 1e-7) {
    throw std::logic_error("Matrix determinant is 0.");
  }
  S21Matrix cofactors = CalcComplements();
  S21Matrix adjugate = cofactors.Transpose();
  double num = 1.0 / det;
  adjugate.MulNumber(num);
  return adjugate;
}

int S21Matrix::get_rows() const { return rows_; }
int S21Matrix::get_cols() const { return cols_; }

void S21Matrix::set_rows(int new_rows) {
  if (new_rows < 0)
    throw std::invalid_argument("Rows count cannot be negative.");
  if (new_rows == rows_) return;
  if (cols_ == 0 && new_rows != 0)
    throw std::invalid_argument("incorect matrix.");

  S21Matrix tmp(new_rows, cols_);

  int min_rows = std::min(rows_, new_rows);
  for (int i = 0; i < min_rows; ++i) {
    std::copy(matrix_[i], matrix_[i] + cols_, tmp.matrix_[i]);
  }

  *this = std::move(tmp);
}

void S21Matrix::set_cols(int new_cols) {
  if (new_cols < 0)
    throw std::invalid_argument("Cols count cannot be negative.");
  if (new_cols == cols_) return;
  if (cols_ == 0 && new_cols != 0)
    throw std::invalid_argument("incorect matrix.");

  S21Matrix tmp(rows_, new_cols);

  int min_rows = std::min(cols_, new_cols);
  for (int i = 0; i < min_rows; ++i) {
    std::copy(matrix_[i], matrix_[i] + rows_, tmp.matrix_[i]);
  }
  *this = std::move(tmp);

}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
  return result;
}


bool S21Matrix::operator==(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(other.matrix_[i][j] - matrix_[i][j]) > 1e-7) {
        return false;
      }
    }
  }
  return true;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }

  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  if (rows_ != other.get_rows() || cols_ != other.get_cols()) {
    throw std::logic_error("different matrix dimensions.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }

  return *this;
}

S21Matrix& S21Matrix::operator*=(const double number) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= number;
    }
  }

  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  if (cols_ != other.get_rows()) {
    throw std::logic_error(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }
  S21Matrix result{rows_, other.get_cols()};
  for (int i = 0; i < result.get_rows(); ++i) {
    for (int j = 0; j < result.get_cols(); ++j) {
      for (int k = 0; k < cols_; ++k) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(result);

  return *this;
}

double& S21Matrix::operator()(int i, int j) { return matrix_[i][j]; }


// не нужный алгоритм Штрассена (потому что я напридумывала хуйни)

int S21Matrix::NextPow2(int x) {
  int p = 1;
  while (p < x) p <<=1;
    return p;
}

S21Matrix S21Matrix::PadToSquare(int n) const {
    S21Matrix res(n, n);
    for (int i = 0; i < rows_; ++i) {
    	for (int j = 0; j < cols_; ++j) {
         res.matrix_[i][j] = matrix_[i][j];
        }
  	}
  	return res;
}

S21Matrix S21Matrix::Slice(int r0, int c0, int sz) const {
  	S21Matrix res(sz, sz);
  	for (int i = 0; i < sz; ++i) {
    	for (int j = 0; j < sz; ++j) {
      	res.matrix_[i][j] = matrix_[r0 + i][c0 + j];
    	}
  	}
 	return res;
}

void S21Matrix::Paste(const S21Matrix& block, int r0, int c0) {
  for (int i = 0; i < block.rows_; ++i) {
    for (int j = 0; j < block.cols_; ++j) {
      matrix_[r0 + i][c0 + j] = block.matrix_[i][j];
    }
  }
}

S21Matrix S21Matrix::StrassenRec(const S21Matrix& A, const S21Matrix& B) {
  int n = A.rows_;
  const int LEAF = 64; // Для оптимизации небольших перемножений (никогда не произойдет)

  if (n <= LEAF) {
    S21Matrix base(A);
    base.MulMatrix(B);
    return base;
  }

  int k = n/2;
  S21Matrix A11 = A.Slice(0, 0, k);
  S21Matrix A12 = A.Slice(0, k, k);
  S21Matrix A21 = A.Slice(k, 0, k);
  S21Matrix A22 = A.Slice(k, k, k);

  S21Matrix B11 = B.Slice(0, 0, k);
  S21Matrix B12 = B.Slice(0, k, k);
  S21Matrix B21 = B.Slice(k, 0, k);
  S21Matrix B22 = B.Slice(k, k, k);

  S21Matrix T1(A11);
  T1.SumMatrix(A22);
  S21Matrix T2(B11);
  T2.SumMatrix(B22);
  S21Matrix M1 = StrassenRec(T1, T2);

  // M2 = (A21 + A22) * B11
  S21Matrix T3(A21);
  T3.SumMatrix(A22);
  S21Matrix M2 = StrassenRec(T3, B11);

  // M3 = A11 * (B12 - B22)
  S21Matrix T4(B12);
  T4.SubMatrix(B22);
  S21Matrix M3 = StrassenRec(A11, T4);

  // M4 = A22 * (B21 - B11)
  S21Matrix T5(B21);
  T5.SubMatrix(B11);
  S21Matrix M4 = StrassenRec(A22, T5);

  // M5 = (A11 + A12) * B22
  S21Matrix T6(A11);
  T6.SumMatrix(A12);
  S21Matrix M5 = StrassenRec(T6, B22);

  // M6 = (A21 - A11) * (B11 + B12)
  S21Matrix T7(A21);
  T7.SubMatrix(A11);
  S21Matrix T8(B11);
  T8.SumMatrix(B12);
  S21Matrix M6 = StrassenRec(T7, T8);

  // M7 = (A12 - A22) * (B21 + B22)
  S21Matrix T9(A12);
  T9.SubMatrix(A22);
  S21Matrix T10(B21);
  T10.SumMatrix(B22);
  S21Matrix M7 = StrassenRec(T9, T10);

  // C11 = M1 + M4 - M5 + M7
  S21Matrix C11(M1);
  C11.SumMatrix(M4);
  C11.SubMatrix(M5);
  C11.SumMatrix(M7);

  // C12 = M3 + M5
  S21Matrix C12(M3);
  C12.SumMatrix(M5);

  // C21 = M2 + M4
  S21Matrix C21(M2);
  C21.SumMatrix(M4);

  // C22 = M1 - M2 + M3 + M6
  S21Matrix C22(M1);
  C22.SubMatrix(M2);
  C22.SumMatrix(M3);
  C22.SumMatrix(M6);

  S21Matrix C(n, n);
  C.Paste(C11, 0, 0);
  C.Paste(C12, 0, k);
  C.Paste(C21, k, 0);
  C.Paste(C22, k, k);

  return C;
}

void S21Matrix::StrassenAlgorithm(const S21Matrix& other) {
  if (cols_ != other.get_rows()) {
    throw std::logic_error(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }

  int m = rows_;
  int n = cols_;
  int p = other.get_cols(); // для конечного результата, определяет, сколько столбцов копировать.

  int s = NextPow2(std::max(m, std::max(n, p)));
  S21Matrix A = PadToSquare(s);
  S21Matrix B = other.PadToSquare(s);

  S21Matrix C = StrassenRec(A, B);
  S21Matrix res(m, p);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < p; ++j) {
      res.matrix_[i][j] = C.matrix_[i][j];
    }
  }
  *this = std::move(res);
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  if (cols_ != other.rows_) {
    throw std::logic_error(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix.");
  }

  S21Matrix result(rows_, other.cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      double sum = 0.0;
      for (int k = 0; k < cols_; ++k) {
        sum += matrix_[i][k] * other.matrix_[k][j];
      }
      result.matrix_[i][j] = sum;
    }
  }

  return result;
}


S21Matrix S21Matrix::operator*(double number) const {
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result[i][j] = matrix_[i][j] * number;
    }
  }
  return result;
}

