#include <complex>
#include <filesystem>
#include <iostream>

#include "material.hpp"
#include "matplotlib.hpp"
#include "matrix.hpp"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <linalg.hpp>

int main()
{
  // Simple Test fmtlib
  fmt::print("Hello World!\n");
  // Test input data
  std::filesystem::path dataFile("../mat/alq3_literature.dat");

  Material alq3 = Material(dataFile, ',');
  std::complex<double> alq3Data = alq3.getRefIndex(4003);
  fmt::print("{}\n", alq3Data.real());
  fmt::print("{}\n", alq3Data.imag());

  // Test slicing. alq3Data contains (wvl, n, k) as | wvl | n | k | wvl | n | k | ...
  // Here we print the first 10 values of wavelength and n
  std::vector<double> x, y;
  linspace(x, -2.0, 2.0, 100);
  for (const auto& val : x) { y.push_back(val * val); }

  // Test Matrix class and type selector specialization
  Matrix<int> m(3, 3);
  for (size_t i = 0; i < m.rows(); ++i) {
    for (size_t j = 0; j < m.cols(); ++j) {
      if (i == j) m(i, j) = 1;
    }
  }

  for (size_t i = 0; i < m.rows(); ++i) {
    for (size_t j = 0; j < m.cols(); ++j) { std::cout << m(i, j) << ", "; }
    std::cout << '\n';
  }

  // Test matplotlib embedding
  figure();
  plot(x, y);
  figure();
  imshow(m);
  // save("./test.png");
  show();
}