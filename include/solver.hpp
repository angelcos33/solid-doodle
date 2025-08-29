#include <complex>
#include <filesystem>
#include <iostream>
#include <vector>
#include <cmath>

#include "material.hpp"
#include "matrix.hpp"
#include "linalg.hpp"



#define PLANCK 6.626e-34
#define PI 3.1415
#define ELEM_CHARGE 1.602e-19
#define RAD_EFF 1
#define TAU0 1e-9 // free space lifetime (seconds)

class MatStack
{
    private:

        //structure-related properties
        unsigned int mDipoleLayer;
        double mLambda;
 
        std::vector<double> mLayerThickness;
        std::vector<Material> mMaterials;

        //descretization 
        long double mXres;
        std::vector<std::complex<long double>> mX;
        std::vector<std::complex<long double>> mU;

        std::vector<std::complex<long double>> mDx;
        std::vector<std::complex<long double>> mDu;

        //wavevector boundaries (might not need it here, I'll check later)
        std::vector<std::complex<double>> mK;

        double mCritU;
        double mRadU;
        double mMatU;

        void setup_discretization();
        void setup_wavevector();
        double big_loop();

    public:
    
        MatStack(unsigned int dipoleLayer, double wavelength, std::vector<double> thicknesses,  std::vector<Material> materials);
};
