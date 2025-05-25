#pragma once

#include <massert.hpp>
#include <utility>
#include <vector>

template<typename T> class Matrix
{

  size_t mRows, mCols;
  std::vector<T> mData;

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

  size_t size() const { return mRows * mCols; }
  size_t rows() const { return mRows; }
  size_t cols() const { return mCols; }
};