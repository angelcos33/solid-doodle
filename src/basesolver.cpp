#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>
#include <Eigen/Core>

#include <forwardDecl.hpp>
#include "basesolver.hpp"
#include "linalg.hpp"

FresnelCoeffs::FresnelCoeffs(const CMatrix& R_perp, const CMatrix& R_para): perp{R_perp}, para{R_para} {}

GFCoeffRatios::GFCoeffRatios(const CMatrix& CB, const CMatrix& FB, const CMatrix& CT, const CMatrix& FT):
                  cb{CB},
                  fb{FB},
                  ct{CT},
                  ft{FT}
                  {}

GFCoeff::GFCoeff(const CMatrix& c1, const CMatrix& c2, const CMatrix& c3, const CMatrix& c4, const CMatrix& c5, const CMatrix& c6):
                 c{c1},
                 cd{c2},
                 f_perp{c3},
                 fd_perp{c4},
                 f_para{c5},
                 fd_para{c6}
                 {}

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
void BaseSolver::calculateFresnelCoeffs(CMatrix& R_perp, CMatrix& R_para) {
    R_perp = (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols()) - matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())) / 
            (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols()) + matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols()));
    (R_perp.bottomRows(R_perp.rows() - mDipoleLayer)) *= -1.0; 

    R_para = ((matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(1, matstack.epsilon.size() - 1) - 
                (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(0, matstack.epsilon.size() - 1));
    R_para /= ((matstack.h.block(0, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(1, matstack.epsilon.size() - 1) + 
                (matstack.h.block(1, 0, matstack.h.rows() - 1, matstack.h.cols())).colwise() * matstack.epsilon.segment(0, matstack.epsilon.size() - 1)); 
    (R_para.bottomRows(R_para.rows() - mDipoleLayer)) *= -1.0; 
}

void BaseSolver::calculateGFCoeffRatios(const FresnelCoeffs& fresnelCoeffs, CMatrix& CB, CMatrix& FB, CMatrix& CT, CMatrix& FT) {
    CMPLX I(0.0, 1.0);
    for (Eigen::Index i = 1; i < mDipoleLayer + 1; ++i) {
        CB.row(i) = Eigen::exp(-2.0 * I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * (fresnelCoeffs.perp.row(i-1) + (CB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1)))) / 
                    (1 + fresnelCoeffs.perp.row(i - 1) * (CB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1))));

        FB.row(i) = Eigen::exp(-2.0 * I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i - 1)) * (-fresnelCoeffs.para.row(i-1) + (FB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1)))) / 
                    (1 - fresnelCoeffs.para.row(i - 1) * (FB.row(i-1) * Eigen::exp(2.0 * I * matstack.h.row(i-1) * (matstack.z0.cast<CMPLX>())(i - 1))));
    }

    for (Eigen::Index i = matstack.numLayers - mDipoleLayer - 2; i >= 0; --i) {
        Eigen::Index indexFromTop = i + mDipoleLayer;
        CT.row(i) = Eigen::exp(2.0 * I * matstack.h.row(indexFromTop) * (matstack.z0.cast<CMPLX>())(indexFromTop)) * (fresnelCoeffs.perp.row(indexFromTop) + (CT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop+1) * (matstack.z0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 + fresnelCoeffs.perp.row(indexFromTop) * (CT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop + 1) * (matstack.z0.cast<CMPLX>())(indexFromTop))));

        FT.row(i) = Eigen::exp(2.0 * I * matstack.h.row(indexFromTop) * (matstack.z0.cast<CMPLX>())(indexFromTop)) * (-fresnelCoeffs.para.row(indexFromTop) + (FT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop+1) * (matstack.z0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 - fresnelCoeffs.para.row(indexFromTop) * (FT.row(i+1) * Eigen::exp(-2.0 * I * matstack.h.row(indexFromTop + 1) * (matstack.z0.cast<CMPLX>())(indexFromTop))));
    }
}

void BaseSolver::calculateGFCoeffs(
                                   const GFCoeffRatios& gfCoeffRatios,
                                   CMatrix& c, 
                                   CMatrix& cd,
                                   CMatrix& f_perp, 
                                   CMatrix& fd_perp, 
                                   CMatrix& f_para, 
                                   CMatrix& fd_para) {

    CMPLX I(0.0, 1.0);
    c.row(mDipoleLayer) = (gfCoeffRatios.cb.row(mDipoleLayer) + 1) * gfCoeffRatios.ct.row(0) / (1 - gfCoeffRatios.cb.row(mDipoleLayer) * gfCoeffRatios.ct.row(0));
    cd.row(mDipoleLayer) = (gfCoeffRatios.ct.row(0) + 1) * gfCoeffRatios.cb.row(mDipoleLayer) / (1 - gfCoeffRatios.cb.row(mDipoleLayer) * gfCoeffRatios.ct.row(0));

    f_perp.row(mDipoleLayer) = (gfCoeffRatios.fb.row(mDipoleLayer) + 1) * gfCoeffRatios.ft.row(0) / (1 - gfCoeffRatios.fb.row(mDipoleLayer) * gfCoeffRatios.ft.row(0));
    fd_perp.row(mDipoleLayer) = (gfCoeffRatios.ft.row(0) + 1) * gfCoeffRatios.fb.row(mDipoleLayer) / (1 - gfCoeffRatios.fb.row(mDipoleLayer) * gfCoeffRatios.ft.row(0));

    f_para.row(mDipoleLayer) = (gfCoeffRatios.fb.row(mDipoleLayer) - 1) * gfCoeffRatios.ft.row(0) / (1 - gfCoeffRatios.fb.row(mDipoleLayer) * gfCoeffRatios.ft.row(0));
    fd_para.row(mDipoleLayer) = (1 - gfCoeffRatios.ft.row(0)) * gfCoeffRatios.fb.row(mDipoleLayer) / (1 - gfCoeffRatios.fb.row(mDipoleLayer) * gfCoeffRatios.ft.row(0));

    Vector boolValue = Vector::Zero(matstack.numLayers);
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
}

void BaseSolver::calculateLifetime(const GFCoeff& gfCoeff, Vector& bPerp, Vector& bPara) {

    CVector bTmp(matstack.u.size() - 1);
    CVector bTmp2(matstack.u.size() - 1);

    bTmp = matstack.dX * (3.0 / 2.0) * Eigen::pow(Eigen::cos(matstack.x.head(matstack.x.size() - 1)), 3);
    bTmp *= (gfCoeff.f_perp(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)) + gfCoeff.fd_perp(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)));
    bPerp = bTmp.real();

    bTmp2 = (gfCoeff.c(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)) + gfCoeff.cd(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)));
    bTmp = Eigen::pow(Eigen::sin(matstack.x.head(matstack.x.size() - 1)), 2); 
    bTmp *= (gfCoeff.f_para(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)) - gfCoeff.fd_para(mDipoleLayer, Eigen::seqN(0, matstack.u.size() - 1)));
    bTmp += bTmp2;
    bTmp *= matstack.dX * (3.0 / 4.0) * Eigen::cos(matstack.x.head(matstack.x.size() - 1));
    bPara = bTmp.real();
}

void BaseSolver::calculateDissPower(const GFCoeff& gfCoeff, const double bPerpSum) {

    // Power calculation
    mPowerPerpU.resize(matstack.numLayers - 1, matstack.u.size() - 1);
    mPowerParaU.resize(matstack.numLayers - 1, matstack.u.size() - 1);
    CMatrix powerParaTemp(matstack.numLayers - 1, matstack.u.size() - 1);

    double q = 1.0; // PLQY
    CMPLX I(0.0, 1.0);

    Vector boolValue = Vector::Zero(matstack.numLayers);
    boolValue(mDipoleLayer) = 1.0;
    for (Eigen::Index i=0; i<matstack.numLayers-1; ++i) {
        mPowerPerpU.row(i) = (-3.0 * q * matstack.dU / 4.0) * ((Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 3)) / Eigen::abs(1 - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) * 
                            (Eigen::sqrt(matstack.epsilon(i) / matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) * 
                            (std::conj(std::sqrt(matstack.epsilon(i))) / std::sqrt(matstack.epsilon(i))); 
        mPowerPerpU.row(i) *= (gfCoeff.f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) - 
                            ((gfCoeff.fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        mPowerPerpU.row(i) *= (Eigen::conj((gfCoeff.f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                            ((gfCoeff.fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)))));         
                            
        mPowerParaU.row(i) = (-3.0 * q * matstack.dU / 8.0) * (matstack.u.segment(0, matstack.u.size() - 1) * Eigen::conj(Eigen::sqrt(matstack.epsilon(i) / matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2)))) /
                            (Eigen::abs(1 - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))); 
        mPowerParaU.row(i) *= (gfCoeff.c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                            ((gfCoeff.cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        mPowerParaU.row(i) *= (Eigen::conj((gfCoeff.c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) -
                            ((gfCoeff.cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)))));
                            
        powerParaTemp.row(i) = (-3.0 * q * matstack.dU / 8.0) * (matstack.u.segment(0, matstack.u.size() - 1) * Eigen::sqrt(matstack.epsilon(i)/matstack.epsilon(mDipoleLayer) - Eigen::pow(matstack.u.segment(0, matstack.u.size() - 1), 2))) *
                            (std::conj(std::sqrt(matstack.epsilon(i))) / std::sqrt(matstack.epsilon(i)));
        powerParaTemp.row(i) *= (gfCoeff.f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) -
                                ((gfCoeff.fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i)));
        powerParaTemp.row(i) *= (Eigen::conj((gfCoeff.f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))) +
                                ((gfCoeff.fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * matstack.h.row(i) * (matstack.z0.cast<CMPLX>())(i))))); 
    }
    mPowerParaU += powerParaTemp;

    // Fraction power calculation
    Matrix m1 = Eigen::real(mPowerPerpU.block(0, 0, mPowerPerpU.rows() - 1, mPowerPerpU.cols()));
    Matrix m2 = Eigen::real(mPowerPerpU.block(1, 0, mPowerPerpU.rows() - 1, mPowerPerpU.cols()));
    mFracPowerPerpU = Eigen::abs(m2 - m1);
    mFracPowerPerpU /= std::abs(bPerpSum);
}

void BaseSolver::calculate() {

    // Loggin
    std::cout << "\n\n\n" << "-----------------------------------------------------------------\n";
    std::cout << "              Starting calculation             \n";
    std::cout << "-----------------------------------------------------------------\n" << "\n\n";

    double q = 1.0; // PLQY. will become a member in the future
    // Matrix and Vector sizes for initialization
    Eigen::Index numInterfaces = matstack.numLayers - 1;
    Eigen::Index numKVectors = matstack.u.size();
    Eigen::Index numLayersTop = mDipoleLayer + 1;
    Eigen::Index numLayersBottom = matstack.numLayers - mDipoleLayer;

    // R_para/perp: Fresnel coefficients
    CMatrix R_perp(numInterfaces, numKVectors);
    CMatrix R_para(numInterfaces, numKVectors);
    // CB and FB coefficients
    CMatrix CB = CMatrix::Zero(numLayersTop, numKVectors);
    CMatrix FB = CMatrix::Zero(numLayersTop, numKVectors);
    // CT and FT coefficients
    CMatrix CT = CMatrix::Zero(numLayersBottom, numKVectors);
    CMatrix FT = CMatrix::Zero(numLayersBottom, numKVectors);
    // c, cd, f, fd coefficients
    CMatrix c = CMatrix::Zero(matstack.numLayers, numKVectors);
    CMatrix cd = CMatrix::Zero(matstack.numLayers, numKVectors);
    CMatrix f_perp = CMatrix::Zero(matstack.numLayers, numKVectors);
    CMatrix fd_perp = CMatrix::Zero(matstack.numLayers, numKVectors);
    CMatrix f_para = CMatrix::Zero(matstack.numLayers, numKVectors);
    CMatrix fd_para = CMatrix::Zero(matstack.numLayers, numKVectors);
    // Dipole lifetime calculations
    Vector bPerp(numKVectors - 1);
    Vector bPara(numKVectors - 1);

    calculateFresnelCoeffs(R_perp, R_para);
    FresnelCoeffs fresnelCoeffs = FresnelCoeffs(R_perp, R_para);

    calculateGFCoeffRatios(fresnelCoeffs, CB, FB, CT, FT);
    GFCoeffRatios gfCoeffRatios = GFCoeffRatios(CB, FB, CT, FT);

    calculateGFCoeffs(gfCoeffRatios, c, cd, f_perp, fd_perp, f_para, fd_para);
    GFCoeff gfCoeff = GFCoeff(c, cd, f_perp, fd_perp, f_para, fd_para);
    
    calculateLifetime(gfCoeff, bPerp, bPara);
    double bPerpSum = 1.0 - q + q * (1 + bPerp.sum());
    double bParaSum = 1.0 - q + q * (1 + bPara.sum());

    calculateDissPower(gfCoeff, bPerpSum);

    // Loggin
    std::cout << "\n\n\n" << "-----------------------------------------------------------------\n";
    std::cout << "              Calculation finished!             \n";
    std::cout << "-----------------------------------------------------------------\n" << "\n\n";
}

void BaseSolver::calculateEmissionSubstrate(Vector& thetaGlass, Vector& powerPerpGlass, Vector& powerParaGlass) {
    double uCriticalGlass = std::real(std::sqrt(matstack.epsilon(matstack.numLayers - 1) / matstack.epsilon(mDipoleLayer)));
    auto uGlassIt = std::find_if(matstack.u.begin(), matstack.u.end(), [uCriticalGlass](auto a){ return a > uCriticalGlass; });
    auto uGlassIndex = uGlassIt - matstack.u.begin();

    thetaGlass = Eigen::real(Eigen::acos(Eigen::sqrt(1 - matstack.epsilon(mDipoleLayer) / matstack.epsilon(matstack.numLayers - 1) * Eigen::pow(matstack.u(Eigen::seq(1, uGlassIndex)), 2))));

    powerPerpGlass = ((Eigen::real(mPowerPerpU(matstack.numLayers-2, Eigen::seq(1, uGlassIndex)))) * std::sqrt(std::real(matstack.epsilon(matstack.numLayers-1) / matstack.epsilon(mDipoleLayer))));
    powerPerpGlass /= Eigen::tan(thetaGlass);

    powerParaGlass = ((Eigen::real(mPowerParaU(matstack.numLayers-2, Eigen::seq(1, uGlassIndex)))) * std::sqrt(std::real(matstack.epsilon(matstack.numLayers-1) / matstack.epsilon(mDipoleLayer))));
    powerParaGlass /= Eigen::tan(thetaGlass);
}

void BaseSolver::modeDissipation(Vector& u, Matrix& fracPowerPerp) {
    u = matstack.u;
    fracPowerPerp = mFracPowerPerpU;
}