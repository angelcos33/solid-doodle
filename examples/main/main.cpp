#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <iterator>
#include <linalg.hpp>
#include <material.hpp>
#include <memory>
#include <numeric>
#include <simulation.hpp>

#include <Eigen/Core>

#include <matplot/matplot.h>

int main()
{
  // Set up stack
  double wavelength = 500;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.emplace_back("/src/mat/air.nk", ',');
  materials.emplace_back("/src/mat/CBP.nk", ',');
  materials.emplace_back("/src/mat/glass.nk", ',');
  materials.emplace_back("/src/mat/glass.nk", ',');

  d.push_back(50e-9);
  d.push_back(5000e-10);

  GaussianSpectrum spectrum(450, 700, wavelength, 1.0);

  // Create Solver
  auto simulation = std::make_unique<SimulationSweep>(materials, d, dipoleLayer, 25e-9, spectrum);

  // Calculate power
  auto start = std::chrono::steady_clock::now();
  simulation->calculate();
  auto finish = std::chrono::steady_clock::now();
  double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
  std::cout << "Elapsed time: " << elapsed_seconds << '\n';

  // Polar figure
  Eigen::ArrayXd thetaGlass, powerPerpAngleGlass, powerParasPolAngleGlass, powerParapPolAngleGlass;
  simulation->calculateEmissionSubstrate(thetaGlass, powerPerpAngleGlass, powerParapPolAngleGlass, powerParasPolAngleGlass);

  std::cout << simulation->mPowerPerpUpPol.leftCols(5) << '\n';

  matplot::figure();
  matplot::plot(thetaGlass, powerPerpAngleGlass, "-o");
  matplot::hold(matplot::on);
  matplot::plot(thetaGlass, powerParapPolAngleGlass, "-o");
  matplot::plot(thetaGlass, powerParasPolAngleGlass, "-o");
  
  matplot::figure();
  matplot::polarplot(thetaGlass, powerPerpAngleGlass, "-")->line_width(2);
  matplot::hold(matplot::on);
  matplot::polarplot(thetaGlass, powerParapPolAngleGlass, "-")->line_width(2);
  matplot::polarplot(thetaGlass, powerParasPolAngleGlass, "-")->line_width(2);

  matplot::show();
}