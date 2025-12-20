#pragma once

#include <vector>
#include <Eigen/Core>

#include "basesolver.hpp"
#include "material.hpp"

class Simulation : public BaseSolver {
    
    void discretize() override;
    
    public:

        Simulation(const std::vector<Material>& materials, const std::vector<double>& thickness, const size_t dipoleLayer, const double dipolePosition, const double wavelength);
        ~Simulation() = default;

        // Make these methods accessible only from Simulation objects. This way we are sured MatStack is properly initialized.
        using BaseSolver::calculateDissPower;
        using BaseSolver::calculateEmissionSubstrate;
        
        //void plot() override;
};