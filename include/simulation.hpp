/*! \file simulation.hpp
    \brief A header file for energy emission simulation.

    The simulation header is used to define the simulation class
    and auxiliar functionalities which are used to simulate the 
    energy emission behavior of the stack in question.
*/
#pragma once

#include <Eigen/Core>
#include <vector>

#include "basesolver.hpp"
#include "material.hpp"

/*! \class Simulation
    \brief A class for energy emission simulation which inherits from baseSolver.

    The Simulation class is used to perform energy emission simulations for the material
    stack of interest. It inherits from the baseSolver base class which acts as an interface 
    for both the Simulation and Fitting classes. Thus, it shares its basic methods and members
    with the Fitting class and includes additional ones especially intented for emissions
    simulation. 
*/
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
      /*!< Simulation class constructor, the constructor takes a (std) vector of class Material containing the materials of the stack to be simulated, 
      a (std) vector of layer thicknesses with matching indices, the index of the dipole layer, the dipole position within the stack and the chosen wavelength 
      to be used for the simulation */

    ~Simulation() = default;

    void calculateEmissionSubstrate(Vector& thetaGlass, Vector& powerPerpGlass, Vector& powerParaGlass);
    /*!< Member method of Simulation used to calculate the emitted power leaving the substrate. It fills the first position argument, an Eigen ArrayXd object containing the angles, 
    with the angles up to the critical angle (of total reflection) as determined by the characteristics of the stack's final layer. The second and third positional arguments
    are filled with the perpendicular and parallel components of the emitted power. */

    // Make these methods accessible only from Simulation objects. This way we are sure MatStack is properly initialized.
    using BaseSolver::calculate;
    using BaseSolver::modeDissipation;

    // void plot() override;
};