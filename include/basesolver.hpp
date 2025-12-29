#pragma once

#include <vector>

#include <Eigen/Core>
#include "material.hpp"
#include <forwardDecl.hpp>

struct FresnelCoeffs {
    const CMatrix& perp;
    const CMatrix& para;
    FresnelCoeffs(const CMatrix& R_perp, const CMatrix& R_para);
};

struct GFCoeffRatios {
    const CMatrix& cb;
    const CMatrix& fb;
    const CMatrix& ct;
    const CMatrix& ft;
    GFCoeffRatios(const CMatrix& CB, const CMatrix& FB, const CMatrix& CT, const CMatrix& FT);
};

struct GFCoeff {
    const CMatrix& c;
    const CMatrix& cd;
    const CMatrix& f_perp;
    const CMatrix& fd_perp;
    const CMatrix& f_para;
    const CMatrix& fd_para;
    GFCoeff(const CMatrix& c1, const CMatrix& c2, const CMatrix& c3, const CMatrix& c4, const CMatrix& c5, const CMatrix& c6);
};

class BaseSolver {
    protected:
        const std::vector<Material>& mMaterials;
        const std::vector<double>& mThickness;
        Eigen::Index mDipoleLayer;
        const double mDipolePosition;
        const double mWvl;
        
        struct MatStack {

            Eigen::Index numLayers;

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

        // Discretization
        virtual void discretize() = 0;

        // Main calculation
        void calculateFresnelCoeffs(CMatrix& R_perp, CMatrix& R_para);
        void calculateGFCoeffRatios(const FresnelCoeffs& fresnelCoeffs, 
                                    CMatrix& CB, 
                                    CMatrix& FB, 
                                    CMatrix& CT, 
                                    CMatrix& FT);
        void calculateGFCoeffs(const GFCoeffRatios& gfCoeffRatios,
                               CMatrix& c, 
                               CMatrix& cd,
                               CMatrix& f_perp, 
                               CMatrix& fd_perp, 
                               CMatrix& f_para, 
                               CMatrix& fd_para);
        void calculateLifetime(const GFCoeff& gfCoeff,
                               Vector& bPerp,
                               Vector& bPara);
        void calculateDissPower(const GFCoeff& gfCoeff, const double bPerpSum);
        void calculate();


        void calculateEmissionSubstrate(Vector& thetaGlass, Vector& powerPerpGlass, Vector& powerParaGlass);
        void modeDissipation(Vector& u, Matrix& fracPowerPerp);
        
        public:
        
            using CMPLX = std::complex<double>;

            BaseSolver(const std::vector<Material>& materials, const std::vector<double>& thickness, const size_t dipoleLayer, const double dipolePosition, const double wavelength);
            virtual ~BaseSolver() = default;

            CMatrix mPowerPerpU;
            CMatrix mPowerParaU;

            Matrix mFracPowerPerpU;

            //virtual void plot()=0;

    };