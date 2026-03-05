#pragma once

#include <Eigen/Core>
#include <vector>

#include "basesolver.hpp"
#include "material.hpp"

class Simulation : public BaseSolver
{

  void loadMaterialData() override;
  void genInPlaneWavevector() override;
  void genOutofPlaneWavevector() override;
  void discretize() override;

  public:
    Simulation(const std::vector<Material>& materials,
      const std::vector<double>& thickness,
      const size_t dipoleLayer,
      const double dipolePosition,
      const double wavelength);
    ~Simulation() = default;

    void calculateEmissionSubstrate(Vector& thetaGlass, Vector& powerPerpGlass, Vector& powerParaGlass);

    // Make these methods accessible only from Simulation objects. This way we are sured MatStack is properly initialized.
    using BaseSolver::calculate;
    using BaseSolver::modeDissipation;

    // void plot() override;
};