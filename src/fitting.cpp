#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "files.hpp"
#include "fitting.hpp"
#include "linalg.hpp"
#include "material.hpp"
#include "simulation.hpp"

void Fitting::loadMaterialData() {
  // Logging
  std::cout << "\n\n\n"
            << "-----------------------------------------------------------------\n";
  std::cout << "              Loading material data(fitting mode)             \n";
  std::cout << "-----------------------------------------------------------------\n"
            << "\n\n";

  matstack.numLayers = static_cast<Eigen::Index>(mMaterials.size());
  matstack.epsilon.resize(matstack.numLayers);
  for (size_t i = 0; i < static_cast<size_t>(matstack.numLayers); ++i) {
    matstack.epsilon(i) = mMaterials[i].getEpsilon(mWvl);
    std::cout << "Layer " << i << "; Material: (" << matstack.epsilon(i).real() << ", " << matstack.epsilon(i).imag()
              << ")\n";
  }
  
  //getting sim data
  size_t counter = 0;
  Vector xData(mIntensityData.size());
  Vector expIntensity(mIntensityData.size());
  for(std::map<double, double>::iterator it = mIntensityData.begin(); it!=mIntensityData.end(); ++it) {
    xData(i) = it->first;
    expIntensity(i) = it->second;
    counter++;
  }
  matstack.x.resize(xData.size());
  matstack.x.head(xData.size()) = xData.cast<CMPLX>();
}

void Fitting::genInPlaneWavevector() {
  // Cumulative sum of thicknesses
  matstack.z0.resize(matstack.numLayers - 1);
  matstack.z0(0) = 0.0;
  std::partial_sum(mThickness.begin(), mThickness.end(), std::next(matstack.z0.begin()), std::plus<double>());
  matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition);

  // Discretization of in-plane wavevector
  matstack.u = matstack.x.cos().real();

  // Differences
  matstack.dU = matstack.u.segment(1, matstack.u.size() - 1) - matstack.u.segment(0, matstack.u.size() - 1);
  matstack.dX = matstack.x.segment(1, matstack.x.size() - 1) - matstack.x.segment(0, matstack.x.size() - 1);
}

void Fitting::genOutofPlaneWavevector() {
  // Out of plane wavevector
  matstack.k = 2 * M_PI / mWvl / 1e-9 * matstack.epsilon.sqrt();
  matstack.h.resize(matstack.numLayers, matstack.u.size());
  matstack.h = matstack.k(mDipoleLayer) *
               (((matstack.epsilon.replicate(1, matstack.x.size())) / matstack.epsilon(mDipoleLayer)).rowwise() -
                 matstack.x.cos().pow(2).transpose())
                 .sqrt();
}

void Fitting::discretize() {
  loadMaterialData();
  genInPlaneWavevector();
  genOutofPlaneWavevector();
}

Fitting::Fitting(const std::vector<Material>& materials,
  const std::vector<double>& thickness,
  const size_t dipoleLayer,
  const double dipolePosition,
  const double wavelength,
  const std::map<double, double>& expData) :
  BaseSolver(materials,
    thickness,
    dipoleLayer,
    dipolePosition,
    wavelength), // initialization must be performed this way due to const members
  mIntensityData{expData}
{
  if (materials.size() != thickness.size() + 2) {
    throw std::runtime_error("Invalid Input! Number of materials different than number of layers.");
  }
  if (dipoleLayer >= materials.size()) { throw std::runtime_error("Invalid Input! Dipole position is out of bounds."); }
  // Log initialization of Simulation
  std::cout << "\n\n\n"
            << "-----------------------------------------------------------------\n";
  std::cout << "              Initializing Fitting             \n";
  std::cout << "-----------------------------------------------------------------\n"
            << "\n\n";
  discretize();
};