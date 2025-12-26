#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>
#include <Eigen/Core>

#include "linalg.hpp"
#include "simulation.hpp"
#include "material.hpp"

void Simulation::loadMaterialData() {
    matstack.numLayers = static_cast<Eigen::Index>(mMaterials.size());
    matstack.epsilon.resize(matstack.numLayers);
    for (size_t i=0; i<static_cast<size_t>(matstack.numLayers); ++i) {
        matstack.epsilon(i) = mMaterials[i].getEpsilon(mWvl);
        std::cout << "Layer " << i << "; Material: (" << matstack.epsilon(i).real() << ", " << matstack.epsilon(i).imag() << ")\n";
    }
}

void Simulation::genInPlaneWavevector() {
    // Cumulative sum of thicknesses
    matstack.z0.resize(matstack.numLayers - 1);
    matstack.z0(0) = 0.0;
    std::partial_sum(mThickness.begin(), mThickness.end(), std::next(matstack.z0.begin()), std::plus<double>());
    matstack.z0 -= (matstack.z0(mDipoleLayer - 1) + mDipolePosition); 

    // Discretization of in-plane wavevector
    CMPLX I(0.0, 1.0);
    double x_res = 5e-4;
    Vector x_real = arange<Vector>(-M_PI_2, -x_res, x_res);
    CVector x_imag = I * arange<Vector>(x_res, 1.8, x_res);
    matstack.x.resize(x_real.rows() + x_imag.rows());
    matstack.x.head(x_real.size()) = x_real.cast<CMPLX>();
    matstack.x.segment(x_real.size(), x_imag.size()) = x_imag;
    matstack.u = matstack.x.cos().real();

    // Differences
    matstack.dU = matstack.u.segment(1, matstack.u.size() - 1) - matstack.u.segment(0, matstack.u.size() - 1);
    matstack.dX = matstack.x.segment(1, matstack.x.size() - 1) - matstack.x.segment(0, matstack.x.size() - 1);
}

void Simulation::genOutofPlaneWavevector() {
    // Out of plane wavevector
    matstack.k = 2 * M_PI / mWvl / 1e-9 * matstack.epsilon.sqrt();
    matstack.h.resize(matstack.numLayers, matstack.u.size());
    matstack.h = matstack.k(mDipoleLayer) * (((matstack.epsilon.replicate(1, matstack.x.size()))/matstack.epsilon(mDipoleLayer)).rowwise() - matstack.x.cos().pow(2).transpose()).sqrt(); 
}

void Simulation::discretize() {
    loadMaterialData();
    genInPlaneWavevector();
    genOutofPlaneWavevector();
}


Simulation::Simulation(const std::vector<Material>& materials, const std::vector<double>& thickness, const size_t dipoleLayer, const double dipolePosition, const double wavelength):
    BaseSolver(materials, thickness, dipoleLayer, dipolePosition, wavelength) //initialization must be performed this way due to const members
    {
        if (materials.size() != thickness.size() + 2) {
            throw std::runtime_error("Invalid Input! Number of materials different than number of layers.");
        }
        if (dipoleLayer >= materials.size()) {
            throw std::runtime_error("Invalid Input! Dipole position is out of bounds.");
        }
        discretize();
    };