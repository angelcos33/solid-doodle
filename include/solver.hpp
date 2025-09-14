#pragma once

#include <vector>

#include <Eigen/Core>
#include "material.hpp"

class Solver {

    const std::vector<Material> mMaterials;
    const std::vector<double> mThickness;

    Eigen::Index mDipoleLayer;
    const double mDipolePosition;
    const double mWvl;
    Eigen::Index mNumLayers;

    Eigen::ArrayXcd mEpsilon;
    Eigen::ArrayXd mZ0;

    Eigen::ArrayXd mU;    
    Eigen::ArrayXd mdU;    

    Eigen::ArrayXcd mK;
    Eigen::ArrayXXcd mH;

    // Discretization
    void discretize();
    
    public:
        Solver(std::vector<Material> materials, std::vector<double> thickness, size_t dipoleLayer, double dipolePosition, double wavelength);

        Eigen::ArrayXXcd mPowerPerpU;
        Eigen::ArrayXXcd mPowerParaU;

        void calculateDissPower();
        void calculateEmissionSubstrate(Eigen::ArrayXd& thetaGlass, Eigen::ArrayXd& powerGlass);



        

};
