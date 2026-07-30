#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <format>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <linalg.hpp>
#include <material.hpp>
#include <memory>
#include <numeric>
#include <simulation.hpp>

#include <Eigen/Core>

#include <matplot/matplot.h>


void saveToText(const std::string& filename, char delimiter, std::initializer_list<Vector> arrays) {
  Eigen::Index end;
  if (std::all_of(arrays.begin(), arrays.end(), [first = arrays.begin()](const Vector& array){return array.size()==first->size();})) {
   end = arrays.begin()->size();
  }
  else {throw std::runtime_error("All arrays must have same size");}

  std::ofstream file(filename);
  if (file.is_open()) {
    for (Eigen::Index i=0; i < end; ++i) {
      for (const auto& array: arrays) {
        file << std::format("{:5}", array(i)) << delimiter;
      }
      file << std::endl;
    }
  }
  else {throw std::runtime_error("Unable to open file");}
}

int main()
{
  // Set up stack
  const double wavelength = 530;
  std::vector<Material> materials;
  std::vector<double> d;
  const size_t dipoleLayer = 2;

  materials.push_back(Material(0.88, 6.4));
  //materials.push_back(Material(1.0, 0.0));
  materials.push_back(Material(1.9, 0.0));
  materials.push_back(Material(1.78, 0.0));
  materials.push_back(Material(1.52, 0.0));
  materials.push_back(Material(1.91, 0.0));
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
  auto simulation = std::make_unique<Simulation>(SimulationMode::AngleSweep, materials, d, dipoleLayer, 10e-9, spectrum, 0.0, 90.0);

  simulation->calculate();
  // Polar figure
  Eigen::ArrayXd thetaGlass, powerPerpAngleGlass, powerParasPolAngleGlass, powerParapPolAngleGlass;
  simulation->calculateEmissionSubstrate(thetaGlass, powerPerpAngleGlass, powerParapPolAngleGlass, powerParasPolAngleGlass);

  // Save to file
  Vector thetaGlassDeg(thetaGlass.size());
  thetaGlassDeg = thetaGlass * 180 / M_PI;
  //saveToText("OLED_angle_dependent.csv", ',', {thetaGlassDeg, powerPerpAngleGlass, powerParapPolAngleGlass, powerParasPolAngleGlass});

  matplot::figure();
  matplot::plot(thetaGlassDeg, powerPerpAngleGlass, "-o");
  matplot::hold(matplot::on);
  matplot::plot(thetaGlassDeg, powerParapPolAngleGlass, "-o");
  matplot::plot(thetaGlassDeg, powerParasPolAngleGlass, "-o");
  
  matplot::figure();
  matplot::polarplot(thetaGlass, powerPerpAngleGlass, "-")->line_width(2);
  matplot::hold(matplot::on);
  matplot::polarplot(thetaGlass, powerParapPolAngleGlass, "-")->line_width(2);
  matplot::polarplot(thetaGlass, powerParasPolAngleGlass, "-")->line_width(2);

  matplot::show();
}