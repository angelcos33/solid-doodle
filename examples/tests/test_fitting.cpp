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

#include <matplot/matplot.h>

int main() {
  // Set up stack
  double wavelength = 500;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.push_back(Material(wavelength, 1.0, 0.0));
  materials.push_back(Material(wavelength, 1.8, 0.0));
  materials.push_back(Material(wavelength, 1.5, 0.0));
  materials.push_back(Material(wavelength, 1.5, 0.0));

  d.push_back(50e-9);
  d.push_back(5000e-10);

  // Fitting filepath
  const std::string targetToFit("/src/examples/data/setfos_simple_emitter_isotropic.txt");
  // Create Solver
  auto solver = std::make_unique<Fitting>(materials, d, dipoleLayer, 25e-9, wavelength, targetToFit);

  // Fit
  auto fitRes = solver->fitEmissionSubstrate();
  std::cout << fitRes.first << std::endl;
}