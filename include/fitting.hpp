#pragma once

#include <matplot/matplot.h>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <vector>
#include <map>

#include "basesolver.hpp"
#include "material.hpp"


// Generic functor
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

class Fitting : public BaseSolver {

  public:
    //public struct so that the numerical diff struct can access it
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

    ~Fitting() = default;

    Eigen::Array2Xd calculateEmissionSubstrate();
    std::pair<Eigen::VectorXd, Eigen::ArrayXd> Fitting::fitEmissionSubstrate();

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

struct ResFunctorNumericalDiff : Eigen::NumericalDiff<Fitting::ResFunctor>{};