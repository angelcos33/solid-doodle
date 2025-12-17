#include <iostream>
#include <linalg.hpp>
#include <material.hpp>
#include <simulation.hpp>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cmath>
#include <iterator>
#include <numeric>
#include <memory>
#include <functional>
#include <fmt/core.h>
#include <fmt/ranges.h>

#include <Eigen/Core>

#include <matplot/matplot.h>

int main()
{  
  // Set up stack
  double wavelength = 500;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.push_back(Material(wavelength, 1.0, 0.0));
  materials.push_back(Material(wavelength, 1.8, 0.0));
  materials.push_back(Material(wavelength, 2.0, 0.0));
  materials.push_back(Material(wavelength, 1.5, 0.0));
  materials.push_back(Material(wavelength, 1.5, 0.0));

  d.push_back(50e-9);
  d.push_back(20e-9);
  d.push_back(5000e-10);

  // Create Solver
  auto simulation = std::make_unique<Simulation>(materials, d, dipoleLayer, 25e-9, wavelength);

  // Calculate power
  auto start = std::chrono::steady_clock::now();
  simulation->calculateDissPower();
  auto finish = std::chrono::steady_clock::now();
  double elapsed_seconds = std::chrono::duration_cast<
                           std::chrono::duration<double>>(finish - start).count();
  std::cout << "Elapsed time: " << elapsed_seconds << '\n';

  // Polar figure
  Eigen::ArrayXd thetaGlass, powerAngleGlass;
  simulation->calculateEmissionSubstrate(thetaGlass, powerAngleGlass);

  std::cout << simulation->mPowerPerpU.leftCols(5) << '\n';

  matplot::plot(thetaGlass, powerAngleGlass, "-o");
  //matplot::save("test.png");
  matplot::show();
}