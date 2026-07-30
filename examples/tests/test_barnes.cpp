#include <iostream>

#include <material.hpp>
#include <simulation.hpp>

#include <Eigen/Core>
#include <matplot/matplot.h>

/*
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
  Vector const& y = simulation->mFracPowerPerpUpPol.row(0).head(u.size());
  Vector const& yParapPol = simulation->mFracPowerParaUpPol.row(0).head(u.size());
  Vector const& yParasPol = simulation->mFracPowerParaUsPol.row(0).head(u.size());
  std::cout << u.size() << ", " << y.size() << std::endl;
  saveToText("barnes_simulation_results.csv", ',', {u, y, yParapPol, yParasPol});

  matplot::semilogy(u, y)->line_width(2).color("red");
  matplot::hold(matplot::on);
  matplot::semilogy(u, yParapPol)->line_width(2).color("blue");
  matplot::semilogy(u, yParasPol)->line_width(2).color("green");
  matplot::show();
}
*/
int main()
{
  // Set up stack
  const double wavelength = 530;
  std::vector<Material> materials;
  std::vector<double> d;
  const size_t dipoleLayer = 2;

  materials.push_back(Material("/src/mat/Al_Cent.csv", ','));
  materials.push_back(Material(1.9, 0.0));
  materials.push_back(Material("/src/mat/CBP.csv", ','));
  materials.push_back(Material("/src/mat/PEDOT_BaytronP_AL4083.csv", ','));
  materials.push_back(Material("/src/mat/ITO.csv", ','));
  materials.push_back(Material(1.52, 0.0));
  materials.push_back(Material(1.52, 0.0));

  d.push_back(50e-9);
  d.push_back(20e-9);
  d.push_back(35e-9);
  d.push_back(150e-9);
  d.push_back(5000e-10);

  // Spectrum
  const double fwhm = 30;
  GaussianSpectrum spectrum(450, 700, wavelength, fwhm/2.355);

  // Create Solver
  auto simulation = std::make_unique<Simulation>(SimulationMode::ModeDissipation, materials, d, dipoleLayer, 10e-9, spectrum, 0.0, 2.0);

  simulation->calculate();

  // Mode dissipation figure
  Vector const& u = simulation->getInPlaneWavevector();
  Vector const& y = simulation->mFracPowerPerpUpPol.row(1).head(u.size());
  Vector const& yParapPol = simulation->mFracPowerParaUpPol.row(1).head(u.size());
  Vector const& yParasPol = simulation->mFracPowerParaUsPol.row(1).head(u.size());
  //std::cout << u.size() << ", " << y.size() << std::endl;
  //saveToText("OLED_simulation_results.csv", ',', {u, y, yParapPol, yParasPol});

  matplot::semilogy(u, y)->line_width(2).color("red");
  matplot::hold(matplot::on);
  matplot::semilogy(u, yParapPol)->line_width(2).color("blue");
  matplot::semilogy(u, yParasPol)->line_width(2).color("green");
  matplot::xlim({0.0, 2.0});
  matplot::show();
}