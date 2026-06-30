#include <cmath>
#include <filesystem>
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
  double wavelength = 500;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.push_back(Material(1.0, 0.0));
  //materials.push_back(Material(wavelength, 1.8, 0.0));
  materials.emplace_back("/src/mat/CBP.nk", ',');
  materials.push_back(Material(1.5, 0.0));
  materials.push_back(Material(1.5, 0.0));

  d.push_back(50e-9);
  d.push_back(5000e-10);

  // Fitting filepath
  const std::string targetToFit("/src/examples/data/setfos_simple_spectrum_isotropic.txt");
  // Spectrum
  GaussianSpectrum spectrum(450, 700, wavelength, 20.0);
  // Create Solver
  auto solver = std::make_unique<Fitting>(targetToFit, materials, d, dipoleLayer, 25e-9, spectrum);
  // Fit
  auto fitRes = solver->fitEmissionSubstrate();
  std::cout << fitRes.first << std::endl;
}