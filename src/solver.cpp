#include "solver.hpp"

#include <algorithm> 
#include <complex>
#include <filesystem>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#include "material.hpp"
#include "matrix.hpp"
#include "linalg.hpp"

void MatStack::setup_discretization(){
    mXres = 5e-4;

    std::complex<long double> imag{0.0, 1.0};
    std::vector<std::complex<long double>> xReal;
    std::vector<std::complex<long double>> xImag;

    arange<std::complex<long double>>(xReal, -PI/2, -mXres, mXres);
    xReal.push_back(0.0); //I use this to later on take diff([x 0]) into account without creating a new vector

    arange<std::complex<long double>>(xImag, mXres*imag, 1.8, mXres); // why 1.8 though?

    std::vector<std::complex<long double>> x(xReal.size() + xImag.size());
    std::vector<std::complex<long double>> u(x.size());

    std::vector<std::complex<long double>> dx(x.size()-1);
    std::vector<std::complex<long double>> du(x.size()-1);

    std::merge(xReal.begin(), xReal.end(), xImag.begin(), xImag.end(), x);
    std::transform(x.begin(), x.end(), u.begin(), std::cosl);

    for(size_t i = 0; i < xReal.size()-1; i++){ //diff() part of the discretization
        dx[i] = x[i+1] - x[i];
        du[i] = u[i+1] - u[i];
    };
    for(size_t i = xReal.size()-1; i < x.size(); i++){
        dx[i] = -x[i+1] + x[i];
        du[i] = u[i+1] - u[i];
    };

    x.erase(x.begin() + xReal.size()-1); //we erase the additional element here
    u.erase(x.begin() + xReal.size()-1);

    mX = x;
    mU = u;

    mDx = dx;
    mDu = du;
}

void MatStack::setup_wavevector(){
    for(size_t i = 0; i < mMaterials.size(); i++) mK[i] = 2*PI/(mLambda*1e-9*sqrt(mMaterials[i].getEpsilon(mLambda)));
}

MatStack::MatStack(unsigned int dipoleLayer, double wavelength, std::vector<double> thicknesses,  std::vector<Material> materials):
    mDipoleLayer{dipoleLayer},
    mLambda{wavelength},
    mLayerThickness{thicknesses},
    mMaterials{materials},
    mK(mMaterials.size())
    {MatStack::setup_discretization();}

double MatStack::big_loop(){ //I'm coding on a plane, so no stackoverflow to help me out    
}

