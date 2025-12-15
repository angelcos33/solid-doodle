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
        
        //void plot() override;
};