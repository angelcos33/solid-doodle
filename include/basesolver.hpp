#pragma once

#include <vector>

#include <Eigen/Core>
#include "material.hpp"

class BaseSolver {
    protected:
        const std::vector<Material>& mMaterials;
        const std::vector<double>& mThickness;
        Eigen::Index mDipoleLayer;
        const double mDipolePosition;
        const double mWvl;
        
        struct MatStack {

            Eigen::Index numLayers;

            Eigen::ArrayXcd epsilon;
            Eigen::ArrayXd z0;

            Eigen::ArrayXd u;    
            Eigen::ArrayXd dU;    

            Eigen::ArrayXcd k;
            Eigen::ArrayXXcd h;
        };

        MatStack matstack;

        // Discretization
        virtual void discretize() = 0;

        void calculateDissPower();
        void calculateEmissionSubstrate(Eigen::ArrayXd& thetaGlass, Eigen::ArrayXd& powerGlass);
        
        public:
        
            using CMPLX = std::complex<double>;

            BaseSolver(const std::vector<Material>& materials, const std::vector<double>& thickness, const size_t dipoleLayer, const double dipolePosition, const double wavelength);
            virtual ~BaseSolver() = default;

            Eigen::ArrayXXcd mPowerPerpU;
            Eigen::ArrayXXcd mPowerParaU;

            //virtual void plot()=0;

    };