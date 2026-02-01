#pragma once

#include <matplot/matplot.h>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include <vector>
#include <map>

#include "basesolver.hpp"
#include "material.hpp"


// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic> struct Functor
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

/*
struct MyFunctor : Functor<double>
{

  Eigen::ArrayXXd samples;
    
  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
  {
    // x here is vector of fitting params
    for (size_t i = 0; i < this->samples.rows(); ++i) {
      fvec(i) = this->samples(i, 1) -
                (x(0) * std::pow(this->samples(i, 0), 2) + x(1) * this->samples(i, 0) + x(2)); // Errors at each sample
    }
    return 0;
  }

  int inputs() const { return 3; } // Number of fitting params
  int values() const { return this->samples.rows(); } // Number of samples to fit
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFunctor>{};

class Fitting : public BaseSolver
{
  std::map<double, double> mIntensityData;
  
  void loadMaterialData();
  void genInPlaneWavevector();
  void genOutofPlaneWavevector();
  void discretize() override;

  public:
    Fitting(const std::vector<Material>& materials,
      const std::vector<double>& thickness,
      const size_t dipoleLayer,
      const double dipolePosition,
      const double wavelength,
      const map<double, double> expData);
    ~Fitting() = default;

    // Make these methods accessible only from Simulation objects. This way we are sured MatStack is properly initialized.
    using BaseSolver::calculate;
    using BaseSolver::calculateEmissionSubstrate;
    using BaseSolver::modeDissipation;

  // void plot() override;
};
*/