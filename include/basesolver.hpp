/*! \file basesolver.hpp
    \brief A header file which contains the essential calculations needed for both simulating and fitting.

    The basesolver file contains the basesolver class and related structs that perform the fundamental computations 
    needed in order to both fit and simulate the behavior of the stack in question. Its main components are the
    BaseSolver class and the coefficient structs used for the calculations performed by BaseSolver and its child
    classes (Simulation and Fitting).
!*/
#pragma once

#include <vector>

#include "material.hpp"
#include <Eigen/Core>
#include <forwardDecl.hpp>

enum class DipoleDistribution {Uniform};

//! A Struct to contain all the Green's Function coefficients. 
struct SolverCoefficients {
  CMatrix _Rperp{};
  CMatrix _Rpara{};
  CMatrix _cb{};
  CMatrix _fb{};
  CMatrix _ct{};
  CMatrix _ft{};
  CMatrix _c{};
  CMatrix _cd{};
  CMatrix _f_perp{};
  CMatrix _fd_perp{};
  CMatrix _f_para{};
  CMatrix _fd_para{};
};

//!  The BaseSolver (virtual) class. 
/*!
  The BaseSolver virtual class performs the basic crucial calculations needed for both simulation and fitting.
  It is constructed from the properties of the stack in question, namely, a vector of (class) Material instances, 
  a vector of layer thicknesses, the index of the dipole layer, and its position in the stack; and from the wavelength 
  of interest.
  This class provides a common interface for interacting with the fitting and simulation classes. It is meant as a virtual
  base class and thus cannot (and should not) be instatiated directly.
*/
class BaseSolver
{
protected:
  const std::vector<Material>& mMaterials;
  const std::vector<double>& mThickness;
  Eigen::Index mDipoleLayer;
  const double mDipolePosition;
  double mWvl;
  
  //! A struct to represent a stack of materials and its discretization.
  /*! The MatStack struct contains essential information about the stack, such as the distinct points of its discretization, 
  the wavector, its components as well as the permittivities of its materials. This information provides the basis for the 
  calculations performed by its containing class, BaseSolver.
  */
  struct MatStack
  {

    Eigen::Index numLayers;
    Eigen::Index numInterfaces;
    Eigen::Index numLayersTop;
    Eigen::Index numLayersBottom;
    Eigen::Index numKVectors;
 

    CVector epsilon;
    Vector z0;

    CVector x;
    CVector dX;

    Vector u;
    Vector dU;

    CVector k;
    CMatrix h;
  };

  MatStack matstack;
  SolverCoefficients coeffs;

  // Discretization
  void loadMaterialData();
  virtual void discretize() = 0;
  virtual void genInPlaneWavevector() = 0; 
  virtual void genOutofPlaneWavevector() = 0;

  // Main calculation
  void calculateFresnelCoeffs();
  /*!< Function to calculate the fresnel coefficients as a function of the materials inputted. */
  void calculateGFCoeffRatios();
  /*!< Function to calculate the ratios between the coefficients of the dyadic Green functions for both left and right travelling eigenfunctions.*/
  void calculateGFCoeffs();
  /*!< Function to calculate the coefficients of the dyadic Green functions from the ratios obtained from calculateGFcoeffRatios.*/
  void calculateLifetime(Vector& bPerp, Vector& bPara);
  /*!< Function to calculate the lifetime of the dipole*/
  void calculateDissPower(const double bPerpSum);
  /*!< Function to calculate dissipated power at the output. The power is decomposed in its parallel and perpendicular components.*/
  void calculate();
  /*!< Function that initializes that properly initializes all coefficients and call the other member functions sequentially, as needed to 
  obtain the base results needed for both Fitting and Simulation. In particular, the power emitted at the output as given by the real part of 
  the Poynting vector's area integral.*/

  void calculateEmissionSubstrate(Vector& thetaGlass, Vector& powerPerpGlass, Vector& powerParapPolGlass, Vector& powerParasPolGlass) const;
public:
  using CMPLX = std::complex<double>;

  BaseSolver(const std::vector<Material>& materials,
    const std::vector<double>& thickness,
    const size_t dipoleLayer,
    const double dipolePosition,
    const double wavelength);
    /*!< BaseSolver class constructor, the constructor takes a (std) vector of class Material containing the materials of the stack to be simulated, 
    a (std) vector of layer thicknesses with matching indices, the index of the dipole layer, the dipole position within the stack and the chosen wavelength 
    to be used for the essential calculations needed for both Simulation and Fitting.*/

  virtual ~BaseSolver() = default;

  Vector const& getInPlaneWavevector() const;
  
  CMatrix mPowerPerpUpPol;
  CMatrix mPowerParaUpPol;
  CMatrix mPowerParaUsPol;

  Matrix mFracPowerPerpU;
};