#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>
#include <Eigen/Core>

#include "basesolver.hpp"
#include "linalg.hpp"


BaseSolver::BaseSolver(const std::vector<Material>& materials, const std::vector<double>& thickness, const size_t dipoleLayer, const double dipolePosition, const double wavelength):
    mMaterials{materials},
    mThickness{thickness},
    mDipoleLayer{static_cast<Eigen::Index>(dipoleLayer)},
    mDipolePosition{dipolePosition},
    mWvl{wavelength}
    {
        if (materials.size() != thickness.size() + 2) {
            throw std::runtime_error("Invalid Input! Number of materials different than number of layers.");
        }
        if (dipoleLayer >= materials.size()) {
            throw std::runtime_error("Invalid Input! Dipole position is out of bounds.");
        }
    };

void BaseSolver::calculateDissPower() {
    double q = 1.0;
    CMPLX I(0.0, 1.0);

    // R_para/perp: Fresnel coefficients
    Eigen::ArrayXXcd R_perp(matstack.numLayers - 1, matstack.u.size()), R_para(matstack.numLayers - 1, matstack.u.size());
    R_perp = (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols()) - matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())) / 
            (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols()) + matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols()));
    (R_perp.bottomRows(R_perp.rows() - mDipoleLayer)) *= -1.0; 

    R_para = ((matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(1, matstack.epsilon.size() - 1) - 
                (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(0, matstack.epsilon.size() - 1));
    R_para /= ((matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(1, matstack.epsilon.size() - 1) + 
                (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(0, matstack.epsilon.size() - 1)); 
    (R_para.bottomRows(R_para.rows() - mDipoleLayer)) *= -1.0;

    // CB and FB coefficients
    Eigen::ArrayXXcd CB = Eigen::ArrayXXcd::Zero(mDipoleLayer + 1, matstack.u.size());
    Eigen::ArrayXXcd FB = Eigen::ArrayXXcd::Zero(mDipoleLayer + 1, matstack.u.size());
    for (Eigen::Index i = 1; i < mDipoleLayer + 1; ++i) {
        CB.row(i) = Eigen::exp(-2.0 * I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * (R_perp.row(i-1) + (CB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1)))) / 
                    (1 + R_perp.row(i - 1) * (CB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1))));

        FB.row(i) = Eigen::exp(-2.0 * I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * (-R_para.row(i-1) + (FB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1)))) / 
                    (1 - R_para.row(i - 1) * (FB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1))));
    }

    // CT and FT coefficients
    Eigen::ArrayXXcd CT = Eigen::ArrayXXcd::Zero(matstack.numLayers - mDipoleLayer, matstack.u.size());
    Eigen::ArrayXXcd FT = Eigen::ArrayXXcd::Zero(matstack.numLayers - mDipoleLayer, matstack.u.size());
    for (Eigen::Index i = matstack.numLayers - mDipoleLayer - 2; i >= 0; --i) {
        Eigen::Index indexFromTop = i + mDipoleLayer;
        CT.row(i) = Eigen::exp(2.0 * I * matstack.h.row(indexFromTop) * (matstack.z0.cast<CMPLX>())(indexFromTop)) * (R_perp.row(indexFromTop) + (CT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop+1) * (matstack.z0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 + R_perp.row(indexFromTop) * (CT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop + 1) * (matstack.z0.cast<CMPLX>())(indexFromTop))));

        FT.row(i) = Eigen::exp(2.0 * I * matstack.h.row(indexFromTop) * (matstack.z0.cast<CMPLX>())(indexFromTop)) * (-R_para.row(indexFromTop) + (FT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop+1) * (matstack.z0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 - R_para.row(indexFromTop) * (FT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop + 1) * (matstack.z0.cast<CMPLX>())(indexFromTop))));
    }

    // c, cd, f, fd coefficients
    Eigen::ArrayXXcd c = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());
    Eigen::ArrayXXcd cd = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());
    Eigen::ArrayXXcd f_perp = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());
    Eigen::ArrayXXcd fd_perp = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());
    Eigen::ArrayXXcd f_para = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());
    Eigen::ArrayXXcd fd_para = Eigen::ArrayXXcd::Zero(matstack.numLayers, matstack.u.size());

    c.row(mDipoleLayer) = (CB.row(mDipoleLayer) + 1) * CT.row(0) / (1 - CB.row(mDipoleLayer) * CT.row(0));
    cd.row(mDipoleLayer) = (CT.row(0) + 1) * CB.row(mDipoleLayer) / (1 - CB.row(mDipoleLayer) * CT.row(0));

    f_perp.row(mDipoleLayer) = (FB.row(mDipoleLayer) + 1) * FT.row(0) / (1 - FB.row(mDipoleLayer) * FT.row(0));
    fd_perp.row(mDipoleLayer) = (FT.row(0) + 1) * FB.row(mDipoleLayer) / (1 - FB.row(mDipoleLayer) * FT.row(0));

    f_para.row(mDipoleLayer) = (FB.row(mDipoleLayer) - 1) * FT.row(0) / (1 - FB.row(mDipoleLayer) * FT.row(0));
    fd_para.row(mDipoleLayer) = (1 - FT.row(0)) * FB.row(mDipoleLayer) / (1 - FB.row(mDipoleLayer) * FT.row(0));

    Eigen::ArrayXd boolValue = Eigen::ArrayXd::Zero(matstack.numLayers);
    boolValue(mDipoleLayer) = 1.0;

    for (Eigen::Index i = mDipoleLayer; i>=1; --i) {
        c.row(i - 1) = 0.5 * Eigen::exp(I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + c.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (1 + matstack.h.row(i) / matstack.h.row(i - 1)) + 
                    cd.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (1 - matstack.h.row(i) / matstack.h.row(i - 1)));

        cd.row(i - 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + c.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (1 - matstack.h.row(i) / matstack.h.row(i - 1)) + 
                    cd.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (1 + matstack.h.row(i) / matstack.h.row(i - 1)));

        f_perp.row(i - 1) = 0.5 * Eigen::exp(I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) + matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)) + 
                            fd_perp.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) - matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)));

        fd_perp.row(i - 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) - matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)) + 
                            fd_perp.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) + matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)));

        f_para.row(i - 1) = 0.5 * Eigen::exp(I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (matstack.k(i) / matstack.k(i - 1) + matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)) + 
                    fd_para.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                    (matstack.k(i) / matstack.k(i - 1) - matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)));

        fd_para.row(i - 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i - 1) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) - matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)) + 
                            fd_para.row(i) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * 
                            (matstack.k(i) / matstack.k(i - 1) + matstack.k(i - 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i - 1)));
    }

    for (Eigen::Index i = mDipoleLayer; i < matstack.numLayers - 1; ++i) {
        c.row(i + 1) = 0.5 * Eigen::exp(I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (c.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (1 + matstack.h.row(i) / matstack.h.row(i + 1)) + 
                    (boolValue(i) + cd.row(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (1 - matstack.h.row(i) / matstack.h.row(i + 1)));

        cd.row(i + 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (c.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (1 - matstack.h.row(i) / matstack.h.row(i + 1)) + 
                    (boolValue(i) + cd.row(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (1 + matstack.h.row(i) / matstack.h.row(i + 1)));

        f_perp.row(i + 1) = 0.5 * Eigen::exp(I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (f_perp.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) + matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)) + 
                            (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) - matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)));

        fd_perp.row(i + 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (f_perp.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) - matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)) + 
                            (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) + matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)));

        f_para.row(i + 1) = 0.5 * Eigen::exp(I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (f_para.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (matstack.k(i) / matstack.k(i + 1) + matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)) + 
                    (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                    (matstack.k(i) / matstack.k(i + 1) - matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)));

        fd_para.row(i + 1) = 0.5 * Eigen::exp(-I * matstack.h.row(i + 1) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (f_para.row(i) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) - matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)) + 
                            (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)) * 
                            (matstack.k(i) / matstack.k(i + 1) + matstack.k(i + 1) / matstack.k(i) * matstack.h.row(i) / matstack.h.row(i + 1)));
    }

    mPowerPerpU.resize(matstack.numLayers - 1, matstack.u.size() - 1);
    mPowerParaU.resize(matstack.numLayers - 1, matstack.u.size() - 1);
    Eigen::ArrayXXcd powerParaTemp(matstack.numLayers - 1, matstack.u.size() - 1);

    for (Eigen::Index i=0; i<matstack.numLayers-1; ++i) {
        mPowerPerpU.row(i) = (-3.0 * q * matstack.dU / 4.0) * ((Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 3)) / Eigen::abs(1 - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) * 
                            (Eigen::sqrt(matstack.epsilon(i) / matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) * 
                            (std::conj(std::sqrt(matstack.epsilon(i))) / std::sqrt(matstack.epsilon(i))); 
        mPowerPerpU.row(i) *= (f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) - 
                            ((fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        mPowerPerpU.row(i) *= (Eigen::conj((f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                            ((fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)))));         
                            
        mPowerParaU.row(i) = (-3.0 * q * matstack.dU / 8.0) * (matstack.u.segment(0, matstack.u.size() - 1) * Eigen::conj(Eigen::sqrt(matstack.epsilon(i) / matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2)))) /
                            (Eigen::abs(1 - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))); 
        mPowerParaU.row(i) *= (c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                            ((cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        mPowerParaU.row(i) *= (Eigen::conj((c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) -
                            ((cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)))));
                            
        powerParaTemp.row(i) = (-3.0 * q * matstack.dU / 8.0) * (matstack.u.segment(0, matstack.u.size() - 1) * Eigen::sqrt(matstack.epsilon(i)/matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) *
                            (std::conj(std::sqrt(matstack.epsilon(i))) / std::sqrt(matstack.epsilon(i)));
        powerParaTemp.row(i) *= (f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) -
                                ((fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        powerParaTemp.row(i) *= (Eigen::conj((f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                                ((fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))))); 
    }
    mPowerParaU += powerParaTemp;
}

void BaseSolver::calculateEmissionSubstrate(Eigen::ArrayXd& thetaGlass, Eigen::ArrayXd& powerGlass) {
    double uCriticalGlass = std::real(std::sqrt(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer)));
    auto uGlassIt = std::find_if(matstack.u.begin(), matstack.u.end(), [uCriticalGlass](auto a){ return a > uCriticalGlass; });
    auto uGlassIndex = uGlassIt - matstack.u.begin();

    thetaGlass = Eigen::real(Eigen::acos(Eigen::sqrt(1 - matstack.epsilon(mDipoleLayer) / matstack.epsilon(matstack.numLayers - 1) * Eigen::pow(matstack.u(Eigen::seq(1, uGlassIndex)), 2))));
    powerGlass = ((Eigen::real(mPowerPerpU(matstack.numLayers-2, Eigen::seq(1, uGlassIndex)))) * std::sqrt(std::real(matstack.epsilon(matstack.numLayers-1) / matstack.epsilon(mDipoleLayer))));
    powerGlass /= Eigen::tan(thetaGlass);
}