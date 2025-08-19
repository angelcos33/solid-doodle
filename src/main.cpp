#include <complex>
#include <filesystem>
#include <iostream>

#include "material.hpp"
#include "matplotlib.hpp"
#include "matrix.hpp"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <complex.h>
#include <linalg.hpp>

int main()
{
  // Simple Test fmtlib
  fmt::print("Hello World!\n");
  // Test input data
  std::filesystem::path dataFile("../mat/alq3_literature.dat");

  Material alq3 = Material(dataFile, ',');
  std::complex<double> alq3Data = alq3.getRefIndex(4003);
  fmt::print("n: {}\n", alq3Data.real());
  fmt::print("k: {}\n", alq3Data.imag());
  alq3Data = alq3.getEpsilon(4003);
  fmt::print("epsilon1: {}\n", alq3Data.real());
  fmt::print("epsilon2: {}\n", alq3Data.imag());

  // Test slicing. alq3Data contains (wvl, n, k) as | wvl | n | k | wvl | n | k | ...
  // Here we print the first 10 values of wavelength and n
  std::vector<double> x, y;
  linspace(x, -2.0, 2.0, 100);
  y = x * x;

  // Test Matrix class and type selector specialization
  Matrix<int> m(3, 3);
  for (size_t i = 0; i < m.rows(); ++i) {
    for (size_t j = 0; j < m.cols(); ++j) {
      if (i == j) m(i, j) = 1;
    }
  }

  // Test matrix multiplication and constructor with initializer lists
  Matrix<int> m2({{1, 2, 3},
                  {4, 5, 6},
                  {7, 8, 9}});
  Matrix<int> m3 = m2 * m2;

  for (size_t i = 0; i < m3.rows(); ++i) {
    for (size_t j = 0; j < m3.cols(); ++j) { std::cout << m3(i, j) << ", "; }
    std::cout << '\n';
  }

  // Test matplotlib embedding
  figure();
  plot(x, y);
  figure();
  imshow(m3);
  // save("./test.png");
  show();
}