/*! \file fitting.hpp
    \brief A header file for fitting energy emission from experimental data.

    The fitting header is used to define the fitting class
    and auxiliar functionalities which are used to fit the 
    energy emission behavior of the stack in question using experimental data.
*/
#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <vector>
#include <map>

#include "basesolver.hpp"
#include "material.hpp"


// Generic functor
//! A templeted struct used to configure the generic properties of the functor used for optimization. 
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic> 
struct Functor
{
  typedef _Scalar Scalar;
  enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() :
    m_inputs(InputsAtCompileTime),
    m_values(ValuesAtCompileTime)
  {}
  Functor(int inputs, int values) :
    m_inputs(inputs),
    m_values(values)
  {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

/*! \class Fitting
    \brief A class that inherits from BaseSolver to fit experimental data.

    The Fitting class is intended to fit energy emission data. It inherits 
    from BaseSolver and thus shares the same base design as the Simulation
    class and includes additional members to implement specific fitting
    functionalities. The fitting is done using the Levenberg-Marquardt algorithm
    from the Eigen library, therefore the fitting process is fast and computationally
    efficient.
*/
class Fitting : public BaseSolver {

  public:
    //public struct so that the numerical diff struct can access it
    /*! \struct ResFunctor
    \brief Struct used as the specific functor needed to pass the objective function to the optimization algorithm.
    */
    struct ResFunctor : Functor<double> {
      Eigen::Array2Xd powerGlass;
      Vector intensities;
      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;

      int inputs() const;
      int outputs() const;
    };

    Fitting(const std::vector<Material>& materials,
      const std::vector<double>& thickness,
      const size_t dipoleLayer,
      const double dipolePosition,
      const double wavelength,
      const std::map<double, double>& expData);
      /*!< Fitting class constructor, the constructor takes a (std) vector of class Material containing the materials of the stack to be simulated, 
      a (std) vector of layer thicknesses with matching indices, the index of the dipole layer, the dipole position within the stack, the chosen wavelength
      and the experimental data to be used for fitting. */

    ~Fitting() = default;

    Eigen::Array2Xd calculateEmissionSubstrate(); //MAKE PRIVATE
    /*!< Member method of Fitting used to simulate the emitted power leaving the substrate. The function simulates parallel and perpendicular components of the power emitted,
    so that it returns an Eigen array where the first and second columnns are the perpendicular and parallel components of the emitted power, repectively. 

    std::pair<Eigen::VectorXd, Eigen::ArrayXd> fitEmissionSubstrate();
    /*!< Member method of Fitting used to fit the emitted power leaving the substrate. The function uses the components of the power emitted simulated by
    calculateEmissionSubstrate() and the experimentally obtained intensities in order to compute the residuals for fitting. It uses the Levenberg-Marquadt algorithm to 
    optimize the fittinng parameters and returns a (std) pair containing an Eigen vector of optimized parameters and the Eigen array of emitted power as a function of angle.*/

    // Make these methods accessible only from Simulation objects. This way we are sured MatStack is properly initialized.
    using BaseSolver::calculate;
    using BaseSolver::modeDissipation;

  // void plot() override;

  private:
    std::map<double, double> mIntensityData;
    Vector mIntensities;
    Vector mThetaData;

    ResFunctor mResidual;

    void loadMaterialData() override;
    void genInPlaneWavevector() override; 
    void genOutofPlaneWavevector() override;
    void discretize() override;

};

/*! \struct ResFunctorNumericalDiff
    \brief Struct used to calculate the jacobian for the optimization algorithm.
*/
struct ResFunctorNumericalDiff : Eigen::NumericalDiff<Fitting::ResFunctor>{};