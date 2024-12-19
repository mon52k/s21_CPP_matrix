#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : S21Matrix(1, 1) {}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument(
        "The dimension of the matrix cannot be negative or equal to zero");
  } else {
    rows_ = rows;
    cols_ = cols;
    matrix_ = new double*[rows_];
    for (int i = 0; i != rows_; ++i) {
      matrix_[i] = new double[cols_];
      for (int j = 0; j != cols_; ++j) {
        matrix_[i][j] = 0.;
      }
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i != rows_; ++i) {
    std::memmove(matrix_[i], other.matrix_[i], sizeof(double) * cols_);
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.cols_ = other.rows_ = 0.;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  if (matrix_) {
    for (int i = 0; i != rows_; ++i) {
      if (matrix_[i]) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
  if (rows_ || cols_) {
    rows_ = cols_ = 0.;
  }
}

int S21Matrix::get_rows_() const { return rows_; }

void S21Matrix::set_rows_(int rows) {
  if (rows <= 0) {
    throw std::invalid_argument(
        "The number of rows of the matrix cannot be a negative value");
  } else {
    rows_ = rows;
  }
}

int S21Matrix::get_cols_() const { return cols_; }

void S21Matrix::set_cols_(int cols) {
  if (cols <= 0) {
    throw std::invalid_argument(
        "The number of cols of the matrix cannot be a negative value");
  } else {
    cols_ = cols;
  }
}

double& S21Matrix::operator()(int fst_idx, int lst_idx) {
  if (fst_idx >= rows_ || fst_idx < 0 || lst_idx >= cols_ || lst_idx < 0) {
    throw std::out_of_range("Matrix index out of range");
  } else {
    return matrix_[fst_idx][lst_idx];
  }
}

const double& S21Matrix::operator()(int fst_idx, int lst_idx) const {
  if (fst_idx >= rows_ || fst_idx < 0 || lst_idx >= cols_ || lst_idx < 0) {
    throw std::out_of_range("Matrix index out of range");
  } else {
    return matrix_[fst_idx][lst_idx];
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i != rows_; ++i) {
    for (int j = 0; j != cols_; ++j) {
      if (matrix_[i][j] != other.matrix_[i][j]) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::invalid_argument("Invalid matrices. Cannot add.");
  } else {
    for (int i = 0; i != rows_; ++i) {
      for (int j = 0; j != cols_; ++j) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    throw std::invalid_argument("Invalid matrices. Cannot subtract.");
  } else {
    for (int i = 0; i != rows_; ++i) {
      for (int j = 0; j != cols_; ++j) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i != rows_; ++i) {
    for (int j = 0; j != cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("Invalid matrices. Cannot multiply.");
  } else {
    S21Matrix result(rows_, other.cols_);
    for (int i = 0; i < rows_; i++)
      for (int j = 0; j < other.cols_; j++)
        for (int k = 0; k < cols_; k++)
          result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];

    *this = result;
  }
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix new_matrix(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      new_matrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return new_matrix;
}

S21Matrix S21Matrix::GetMinor(int row, int column, const S21Matrix& A) {
  if (row <= 0 || column <= 0 || A.get_rows_() != A.get_cols_()) {
    throw std::runtime_error("Invalid matrix. Cannot calculate minor.");
  }
  S21Matrix result(A.get_rows_() - 1, A.get_cols_() - 1);
  int tmp_row = 0;
  for (int i = 0; i < A.get_rows_(); ++i) {
    if (i == row - 1) {
      continue;
    }
    int tmp_col = 0;
    for (int j = 0; j < A.get_cols_(); ++j) {
      if (j == column - 1) {
        continue;
      }
      result(tmp_row, tmp_col) = A(i, j);
      ++tmp_col;
    }
    ++tmp_row;
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::runtime_error("Non-square matrix. Cannot compute complements.");
  } else {
    S21Matrix result(rows_, cols_);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        S21Matrix minor = GetMinor(i + 1, j + 1, *this);
        double det = minor.Determinant();
        result.matrix_[i][j] = pow(-1, i + j) * det;
      }
    }

    return result;
  }
}

double S21Matrix::Determinant() {
  if (rows_ != cols_)
    throw std::invalid_argument(
        "Non-square matrix. Cannot compute determinant.");
  if (rows_ == 1)
    return matrix_[0][0];
  else if (rows_ == 2)
    return matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  else {
    double det = 0.0;
    for (int i = 0; i < rows_; ++i) {
      S21Matrix minor = GetMinor(1, i + 1, *this);
      double minor_det = minor.Determinant();
      det += pow(-1, i) * matrix_[0][i] * minor_det;
    }
    return det;
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (rows_ != cols_ || det == 0)
    throw std::runtime_error(
        "Matrix must be square and have a non-zero determinant for inversion.");
  S21Matrix complements = CalcComplements();
  S21Matrix transposed = complements.Transpose();
  S21Matrix inverse = transposed * (1.0 / det);

  return inverse;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  try {
    S21Matrix result(*this);
    result.SumMatrix(other);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  try {
    S21Matrix result(*this);
    result.SubMatrix(other);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  try {
    S21Matrix result = *this;
    result.MulMatrix(other);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix S21Matrix::operator*(const double num) const {
  try {
    S21Matrix result = *this;
    result.MulNumber(num);
    return result;
  } catch (const std::exception& e) {
    throw e;
  }
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  try {
    return this->EqMatrix(other);
  } catch (const std::exception& e) {
    throw e;
  }
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) {
    return *this;
  }
  if (matrix_) {
    for (int i = 0; i != rows_; ++i) {
      if (matrix_[i]) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = new double*[rows_];
  for (int i = 0; i != rows_; ++i) {
    matrix_[i] = new double[cols_];
    memmove(matrix_[i], other.matrix_[i], sizeof(double) * cols_);
  }

  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this == &other) {
    return *this;
  }
  if (matrix_) {
    for (int i = 0; i != rows_; ++i) {
      if (matrix_[i]) {
        delete[] matrix_[i];
      }
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }

  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = other.cols_ = 0.;
  other.matrix_ = nullptr;

  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  try {
    *this = *this + other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  try {
    *this = *this - other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  try {
    *this = *this * other;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  try {
    *this = *this * num;
  } catch (const std::exception& e) {
    throw e;
  }
  return *this;
}