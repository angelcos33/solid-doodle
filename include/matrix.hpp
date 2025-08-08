#pragma once

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <massert.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename T> class Matrix
{

  size_t mRows, mCols;
  std::vector<T> mData;

  void increase_by(const T scalar);
  void increase_by(const Matrix& other);

  void multiply_by(const T scalar);

public:
  Matrix() = default;

  // Initializes a matrix with all 0 entries
  Matrix(size_t nRows, size_t nCols) :
    mRows{nRows},
    mCols{nCols},
    mData(nRows * nCols)
  {}

  Matrix(std::initializer_list<std::initializer_list<T>> iList) {
    m_assert(std::all_of(iList.begin(),
                         iList.end(),
                         [ncols = (iList.begin())->size()](auto l){return l.size() == ncols;}), "All rows must have same number of columns");
    mRows = iList.size();
    mCols = (iList.begin())->size();
    for (const auto& rows: iList) {
      mData.insert(mData.end(), rows.begin(), rows.end());
    }
  }

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

  // Matrix addition overloads
  friend Matrix operator+(Matrix lhs, const Matrix& rhs) {
    lhs += rhs;
    return std::move(lhs);
  }

  friend Matrix& operator+=(Matrix& a, const Matrix& b) {
    a.increase_by(b);
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
void Matrix<T>::increase_by(const Matrix<T>& other) {
  if (this->mRows == other.mRows && this->mCols == other.mCols) {
    for (size_t i = 0; i < mData.size(); ++i) {
      mData[i] += other.mData[i];
    }
  }
  else throw std::runtime_error("Matrix dimension mismatch");
}

template<typename T>
void Matrix<T>::multiply_by(const T scalar) {
  for (auto& val: mData) {
    val *= scalar;
  }
}