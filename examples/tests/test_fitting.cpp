#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <material.hpp>
#include <data.hpp>
#include <fitting.hpp>
#include <simulation.hpp>

#include <matplot/matplot.h>

int main() {
  // Set up stack
  double wavelength = 456;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.push_back(Material(1.0, 0.0));
  //materials.push_back(Material(wavelength, 1.8, 0.0));
  materials.emplace_back(1.7, 0.0);
  materials.push_back(Material(1.5, 0.0));
  materials.push_back(Material(1.5, 0.0));

  d.push_back(35e-9);
  d.push_back(5000e-10);

  // Fitting filepath
  //const std::string targetToFit("/src/examples/data/setfos_simple_spectrum_isotropic.txt");
  const std::string targetToFit("/src/examples/data/3ML_processed.txt");
  // Spectrum
  double fwhm = 20;
  GaussianSpectrum spectrum(400, 700, wavelength, fwhm/2.355);
  // Dipole distribution
  DipoleDistribution dipoleDist(0.0, 35e-9, DipoleDistributionType::Uniform);
  // Create Solver
  auto solver = std::make_unique<Fitting>(targetToFit, materials, d, dipoleLayer, 17.5e-9, spectrum, 0.0, 90.0);
  // Fit
  auto fitRes = solver->fitEmissionSubstrate();
  std::cout << fitRes.first << std::endl;
  
  auto sim = std::make_unique<Simulation>(SimulationMode::AngleSweep, materials, d, dipoleLayer, 17.5e-9, spectrum, 0.0, 90.0);
  sim->calculate();

  Eigen::ArrayXd thetaGlass, powerPerpAngleGlass, powerParasPolAngleGlass, powerParapPolAngleGlass;
  sim->calculateEmissionSubstrate(thetaGlass, powerPerpAngleGlass, powerParapPolAngleGlass, powerParasPolAngleGlass);

  Eigen::ArrayXd res(thetaGlass.size()), angleDeg(thetaGlass.size());
  double scalingFactor = 0.00833604;
  double thetaIP = 0.178204;
  res = scalingFactor * (thetaIP * powerPerpAngleGlass + (1 - thetaIP) * powerParapPolAngleGlass);
  angleDeg = thetaGlass * 180 / M_PI;
  //saveToText("3ML_Fit.txt", angleDeg, res, '\t');


  matplot::figure();
  matplot::plot(thetaGlass, res)->line_width(2);
  matplot::show();
}