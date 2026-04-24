#include <iostream>

#include <material.hpp>
#include <simulation.hpp>

#include <Eigen/Core>
#include <matplot/matplot.h>

int main()
{
  // Set up stack
  double wavelength = 550;
  std::vector<Material> materials;
  std::vector<double> d;
  size_t dipoleLayer = 1;

  materials.push_back(Material(wavelength, 1.0, 0.0));
  materials.push_back(Material(wavelength, 1.578, 0.0));
  materials.push_back(Material(wavelength, 0.0715, 4.1958));
  materials.push_back(Material(wavelength, 0.0715, 4.1958));

  d.push_back(141e-9);
  d.push_back(5000e-10);

  // Create Solver
  auto simulation = std::make_unique<Simulation>(materials, d, dipoleLayer, 0.0, wavelength);

  simulation->calculate();

  // Mode dissipation figure
  Vector const& u = simulation->getInPlaneWavevector();
  Vector const& y = simulation->mFracPowerPerpU.row(0);

  matplot::semilogy(u, y)->line_width(2).color("red");
  matplot::show();
}