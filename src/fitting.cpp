#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <matplot/matplot.h>

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
  Vector thetaData(mIntensityData.size());
  Vector mIntensities(mIntensityData.size());
  for(std::map<double, double>::iterator it = mIntensityData.begin(); it!=mIntensityData.end(); ++it) {
    mThetaData(counter) = it->first;
    mIntensities(counter) = it->second;
    counter++;
}
//QUADRUPLE CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  matstack.u = Eigen::real(Eigen::sqrt(matstack.epsilon(matstack.numLayers - 1)/matstack.epsilon(mDipoleLayer)*(1- pow(Eigen::cos(mThetaData), 2))));
  matstack.x.resize(matstack.u.size());
  matstack.x = matstack.u.acos();
}

void Fitting::genInPlaneWavevector() {
  // Cumulative sum of thicknesses
  matstack.z0.resize(matstack.numLayers - 1);
  matstack.z0(0) = 0.0;
  std::partial_sum(mThickness.begin(), mThickness.end(), std::next(matstack.z0.begin()), std::plus<double>());
  matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition);

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
                 matstack.u.pow(2).transpose()).sqrt();
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
  //setting up functor for fitting
  mResidual.intensities = mIntensities;
  mResidual.powerGlass = calculateEmissionSubstrate();
  }

Eigen::Array2Xd Fitting::calculateEmissionSubstrate() {
  Vector powerPerpGlass;
  Vector powerParaGlass;

  double uCriticalGlass =
    std::real(std::sqrt(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer)));
  auto uGlassIt =
    std::find_if(matstack.u.begin(), matstack.u.end(), [uCriticalGlass](auto a) { return a > uCriticalGlass; });
  auto uGlassIndex = uGlassIt - matstack.u.begin();

  powerPerpGlass = ((Eigen::real(mPowerPerpU(matstack.numLayers - 2, Eigen::seq(1, uGlassIndex)))) *
                    std::sqrt(std::real(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer))));
  powerPerpGlass /= Eigen::tan(mThetaData);

  powerParaGlass = ((Eigen::real(mPowerParaU(matstack.numLayers - 2, Eigen::seq(1, uGlassIndex)))) *
                    std::sqrt(std::real(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer))));
  powerParaGlass /= Eigen::tan(mThetaData);

  Eigen::Array2Xd powerGlass;
  powerGlass.row(0) = powerPerpGlass;
  powerGlass.row(1) = powerParaGlass;
  return powerGlass;
}

int Fitting::ResFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
  // x here is vector of fitting params
  for (size_t i = 0; i < intensities.rows(); ++i) {
    fvec(i) = intensities(i) - (x(0)*powerGlass(0, i) + (1 - x(0))*powerGlass(1, i)); // residual of each sample
  }
  return 0;
} 

int Fitting::ResFunctor::inputs() const {return 1;}

int Fitting::ResFunctor::outputs() const {return this->intensities.rows();}

std::pair<Eigen::VectorXd, Eigen::ArrayXd> Fitting::fitEmissionSubstrate() {
  //returns the vector of parameters and the fitted intensities as a std::pair

  std::vector<double> theta(matstack.x.rows()), yFit(matstack.x.rows()), yExp(mResidual.intensities.rows());

  Eigen::ArrayXd::Map(&theta[0], matstack.x.rows()) = matstack.x.real();
  Eigen::ArrayXd::Map(&yExp[0], mResidual.intensities.rows()) = mResidual.intensities;

  // Setup
  Eigen::VectorXd x(1);
  // Initial guess
  x.fill(0.0);

  ResFunctorNumericalDiff functor;
  Eigen::LevenbergMarquardt<ResFunctorNumericalDiff> lm(functor);
  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1e-10;

  int status = lm.minimize(x);
  std::cout << "Number of iterations: " << lm.iter << '\n';
  std::cout << "Status: " << status << '\n';
  std::cout << "Fitting result: " << x << '\n' << '\n';

  // simulation results
  Eigen::ArrayXd optIntensities(matstack.x.rows());
  optIntensities = x(0)*mResidual.powerGlass.row(0) + (1 - x(0))*mResidual.powerGlass.row(1);
  Eigen::ArrayXd::Map(&yFit[0], matstack.x.rows()) = optIntensities;

  //plotting the results
  matplot::figure();
  
  matplot::scatter(theta, yExp);
  matplot::hold(matplot::on);
  matplot::plot(theta, yFit)->line_width(2).color("red");
  matplot::show();

  return std::pair<Eigen::VectorXd, Eigen::ArrayXd>(x, optIntensities);
};