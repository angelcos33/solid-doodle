#pragma once

#include <massert.hpp>
#include <iostream>
#include <utility>
#include <vector>

template<typename T> class Matrix
{

  size_t mRows, mCols;
  std::vector<T> mData;

  void increase_by(const T scalar);
  void multiply_by(const T scalar);

public:
  Matrix() = default;

  // Initializes a matrix with all 0 entries
  Matrix(size_t nRows, size_t nCols) :
    mRows{nRows},
    mCols{nCols},
    mData(nRows * nCols)
  {}

  T& operator()(size_t row, size_t col)
  {
    m_assert(row < mRows && col < mCols, "Index of out of bound.");
    return mData[col * mCols + row];
  }

  const T& operator()(size_t row, size_t col) const
  {
    m_assert(row < mRows && col < mCols, "Index of out of bound.");
    return mData[col * mCols + row];
  }

  // Scalar addition overloads
  friend Matrix operator+(Matrix a, const T scalar) {
    a += scalar;
    return std::move(a); // To remind that RVO cannot be triggered for function argument a (https://stackoverflow.com/questions/35853204/c-overloading-operator-in-a-template-class)
  }

  friend Matrix operator+(const T scalar, Matrix a) {
    a += scalar;
    return std::move(a);
  }

  friend Matrix& operator+=(Matrix& a, const T scalar) {
    a.increase_by(scalar);
    return a;
  }

  // Scalar multiplication overloads
  friend Matrix operator*(Matrix a, const T scalar) {
    a *= scalar;
    return std::move(a);
  }

  friend Matrix operator*(const T scalar, Matrix a) {
    a *= scalar;
    return std::move(a);
  }

  friend Matrix& operator*=(Matrix& a, const T scalar) {
    a.multiply_by(scalar);
    return a;
  }

  T* data() { return mData.data(); }

  const T* data() const { return mData.data(); }

  size_t size() const { return mRows * mCols; }
  size_t rows() const { return mRows; }
  size_t cols() const { return mCols; }
};

template<typename T>
void Matrix<T>::increase_by(const T scalar) {
  for (auto& val: mData) {
    val += scalar;
  }
}

template<typename T>
void Matrix<T>::multiply_by(const T scalar) {
  for (auto& val: mData) {
    val *= scalar;
  }
}