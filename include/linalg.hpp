/*
Generic library for vector utilities. It is based on the std::vector container.
*/
#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include <massert.hpp>

// Vector arithmetic
// Scalar Multiplication
template<typename T, typename U> std::vector<T> operator*(const std::vector<T>& vec, const U scalar)
{
  std::vector<T> res;
  res.reserve(vec.size());

  for (auto& elem : vec) { res.push_back(scalar * elem); }
  return res;
}

template<typename T, typename U> std::vector<T> operator*(const U scalar, const std::vector<T>& vec)
{
  std::vector<T> res;
  res.reserve(vec.size());

  for (auto& elem : vec) { res.push_back(scalar * elem); }
  return res;
}

// Vector addition
template<typename T> std::vector<T> operator+(const std::vector<T>& lvec, const std::vector<T>& rvec)
{
  m_assert(lvec.size() == rvec.size(), "Only vectors with same size can be added!");
  std::vector<T> res;
  res.reserve(lvec.size());

  for (size_t i = 0; i < lvec.size(); ++i) { res.push_back(lvec[i] + rvec[i]); }
  return res;
}

// Vector element-wise multiplication
template<typename T> std::vector<T> operator*(const std::vector<T>& lvec, const std::vector<T>& rvec)
{
  m_assert(lvec.size() == rvec.size(), "Only vector of same size can be multiplied together element-wise!");
  std::vector<T> res;
  res.reserve(lvec.size());

  for (size_t i = 0; i < lvec.size(); ++i) { res.push_back(lvec[i] * rvec[i]); }
  return res;
}

// Vector slicing
template<typename T> std::vector<T> slice(const std::vector<T>& vec, size_t start, int end, size_t stride = 0)
{
  m_assert((size_t)std::abs(end) <= vec.size(), "Slice upper bound is out of range!");

  std::vector<T> res;
  size_t upperBound = end >= 0 ? end : vec.size() + end;
  res.reserve((upperBound - start) / stride);

  for (size_t i = start; i < upperBound; i += stride) {
    if (i >= vec.size()) break;
    res.push_back(vec[i]);
  }
  return res;
}

template<typename T> void arange(std::vector<T>& out, T start, T stop, T step)
{

  m_assert(stop > start, "Upper bound has to be greater than lower bound!");
  size_t N = ((stop - start) / step) + 1;
  out.reserve(N);

  for (size_t i = 0; i < N; ++i) { out.push_back(start + i * step); }
}

template<typename T> void linspace(std::vector<T>& out, T start, T stop, size_t N)
{
  T step = (stop - start) / (N - 1);
  arange(out, start, stop, step);
}