#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

#include <Eigen/Core>
#include "linalg.hpp"
#include "solver.hpp"

using CMPLX = std::complex<double>;

Solver::Solver(std::vector<Material> materials, std::vector<double> thickness, size_t dipoleLayer, double dipolePosition, double wavelength):
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
        Solver::discretize();
    }

void Solver::discretize() {
    mNumLayers = static_cast<Eigen::Index>(mMaterials.size());
    mEpsilon.resize(mNumLayers);
    for (size_t i=0; i<static_cast<size_t>(mNumLayers); ++i) {
        mEpsilon(i) = mMaterials[i].getEpsilon(mWvl);
    }


    // Cumulative sum
    mZ0.resize(mNumLayers - 1);
    mZ0(0) = 0.0;
    std::partial_sum(mThickness.begin(), mThickness.end(), std::next(mZ0.begin()), std::plus<double>());
    mZ0 -= (mZ0(mDipoleLayer - 1) + mDipolePosition);

    // Discretization of in-plane wavevector
    CMPLX I(0.0, 1.0);
    double x_res = 5e-4;
    Eigen::ArrayXd x_real = arange<Eigen::ArrayXd>(-M_PI_2, -x_res, x_res);
    Eigen::ArrayXcd x_imag = I * arange<Eigen::ArrayXd>(x_res, 1.8, x_res);
    Eigen::ArrayXcd x(x_real.rows() + x_imag.rows());
    x.head(x_real.size()) = x_real.cast<CMPLX>();
    x.segment(x_real.size(), x_imag.size()) = x_imag;
    mU = x.cos().real();

    // Differences
    mdU = mU.segment(1, mU.size() - 1) - mU.segment(0, mU.size() - 1);

    // Out of plane wavevector
    mK = 2 * M_PI / mWvl / 1e-9 * mEpsilon.sqrt();
    mH.resize(mNumLayers, x.size());
    mH = mK(mDipoleLayer) * (((mEpsilon.replicate(1, x.size()))/mEpsilon(mDipoleLayer)).rowwise() - x.cos().pow(2).transpose()).sqrt();
}

void Solver::calculateDissPower() {
    double q = 1.0;
    CMPLX I(0.0, 1.0);

    // R_para/perp: Fresnel coefficients
    Eigen::ArrayXXcd R_perp(mNumLayers - 1, mU.size()), R_para(mNumLayers - 1, mU.size());
    R_perp = (mH.block(1, 0, mH.rows() - 1, mH.cols()) - mH.block(0, 0, mH.rows() - 1, mH.cols())) / 
            (mH.block(1, 0, mH.rows() - 1, mH.cols()) + mH.block(0, 0, mH.rows() - 1, mH.cols()));
    (R_perp.bottomRows(R_perp.rows() - mDipoleLayer)) *= -1.0; 

    R_para = ((mH.block(0, 0, mH.rows() - 1, mH.cols())).colwise() * mEpsilon.segment(1, mEpsilon.size() - 1) - 
                (mH.block(1, 0, mH.rows() - 1, mH.cols())).colwise() * mEpsilon.segment(0, mEpsilon.size() - 1));
    R_para /= ((mH.block(0, 0, mH.rows() - 1, mH.cols())).colwise() * mEpsilon.segment(1, mEpsilon.size() - 1) + 
                (mH.block(1, 0, mH.rows() - 1, mH.cols())).colwise() * mEpsilon.segment(0, mEpsilon.size() - 1)); 
    (R_para.bottomRows(R_para.rows() - mDipoleLayer)) *= -1.0;

    // CB and FB coefficients
    Eigen::ArrayXXcd CB = Eigen::ArrayXXcd::Zero(mDipoleLayer + 1, mU.size());
    Eigen::ArrayXXcd FB = Eigen::ArrayXXcd::Zero(mDipoleLayer + 1, mU.size());
    for (Eigen::Index i = 1; i < mDipoleLayer + 1; ++i) {
        CB.row(i) = Eigen::exp(-2.0 * I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * (R_perp.row(i-1) + (CB.row(i-1) * Eigen::exp(2.0 * I * mH.row(i-1) * (mZ0.cast<CMPLX>())(i - 1)))) / 
                    (1 + R_perp.row(i - 1) * (CB.row(i-1) * Eigen::exp(2.0 * I * mH.row(i-1) * (mZ0.cast<CMPLX>())(i - 1))));

        FB.row(i) = Eigen::exp(-2.0 * I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * (-R_para.row(i-1) + (FB.row(i-1) * Eigen::exp(2.0 * I * mH.row(i-1) * (mZ0.cast<CMPLX>())(i - 1)))) / 
                    (1 - R_para.row(i - 1) * (FB.row(i-1) * Eigen::exp(2.0 * I * mH.row(i-1) * (mZ0.cast<CMPLX>())(i - 1))));
    }

    // CT and FT coefficients
    Eigen::ArrayXXcd CT = Eigen::ArrayXXcd::Zero(mNumLayers - mDipoleLayer, mU.size());
    Eigen::ArrayXXcd FT = Eigen::ArrayXXcd::Zero(mNumLayers - mDipoleLayer, mU.size());
    for (Eigen::Index i = mNumLayers - mDipoleLayer - 2; i >= 0; --i) {
        Eigen::Index indexFromTop = i + mDipoleLayer;
        CT.row(i) = Eigen::exp(2.0 * I * mH.row(indexFromTop) * (mZ0.cast<CMPLX>())(indexFromTop)) * (R_perp.row(indexFromTop) + (CT.row(i+1) * Eigen::exp(-2.0 * I * mH.row(indexFromTop+1) * (mZ0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 + R_perp.row(indexFromTop) * (CT.row(i+1) * Eigen::exp(-2.0 * I * mH.row(indexFromTop + 1) * (mZ0.cast<CMPLX>())(indexFromTop))));

        FT.row(i) = Eigen::exp(2.0 * I * mH.row(indexFromTop) * (mZ0.cast<CMPLX>())(indexFromTop)) * (-R_para.row(indexFromTop) + (FT.row(i+1) * Eigen::exp(-2.0 * I * mH.row(indexFromTop+1) * (mZ0.cast<CMPLX>())(indexFromTop)))) / 
                    (1 - R_para.row(indexFromTop) * (FT.row(i+1) * Eigen::exp(-2.0 * I * mH.row(indexFromTop + 1) * (mZ0.cast<CMPLX>())(indexFromTop))));
    }

    // c, cd, f, fd coefficients
    Eigen::ArrayXXcd c = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());
    Eigen::ArrayXXcd cd = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());
    Eigen::ArrayXXcd f_perp = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());
    Eigen::ArrayXXcd fd_perp = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());
    Eigen::ArrayXXcd f_para = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());
    Eigen::ArrayXXcd fd_para = Eigen::ArrayXXcd::Zero(mNumLayers, mU.size());

    c.row(mDipoleLayer) = (CB.row(mDipoleLayer) + 1) * CT.row(0) / (1 - CB.row(mDipoleLayer) * CT.row(0));
    cd.row(mDipoleLayer) = (CT.row(0) + 1) * CB.row(mDipoleLayer) / (1 - CB.row(mDipoleLayer) * CT.row(0));

    f_perp.row(mDipoleLayer) = (FB.row(mDipoleLayer) + 1) * FT.row(0) / (1 - FB.row(mDipoleLayer) * FT.row(0));
    fd_perp.row(mDipoleLayer) = (FT.row(0) + 1) * FB.row(mDipoleLayer) / (1 - FB.row(mDipoleLayer) * FT.row(0));

    f_para.row(mDipoleLayer) = (FB.row(mDipoleLayer) - 1) * FT.row(0) / (1 - FB.row(mDipoleLayer) * FT.row(0));
    fd_para.row(mDipoleLayer) = (1 - FT.row(0)) * FB.row(mDipoleLayer) / (1 - FB.row(mDipoleLayer) * FT.row(0));

    Eigen::ArrayXd boolValue = Eigen::ArrayXd::Zero(mNumLayers);
    boolValue(mDipoleLayer) = 1.0;

    for (Eigen::Index i = mDipoleLayer; i>=1; --i) {
        c.row(i - 1) = 0.5 * Eigen::exp(I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + c.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (1 + mH.row(i) / mH.row(i - 1)) + 
                    cd.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (1 - mH.row(i) / mH.row(i - 1)));

        cd.row(i - 1) = 0.5 * Eigen::exp(-I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + c.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (1 - mH.row(i) / mH.row(i - 1)) + 
                    cd.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (1 + mH.row(i) / mH.row(i - 1)));

        f_perp.row(i - 1) = 0.5 * Eigen::exp(I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) + mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)) + 
                            fd_perp.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) - mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)));

        fd_perp.row(i - 1) = 0.5 * Eigen::exp(-I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) - mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)) + 
                            fd_perp.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) + mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)));

        f_para.row(i - 1) = 0.5 * Eigen::exp(I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (mK(i) / mK(i - 1) + mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)) + 
                    fd_para.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                    (mK(i) / mK(i - 1) - mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)));

        fd_para.row(i - 1) = 0.5 * Eigen::exp(-I * mH.row(i - 1) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) - mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)) + 
                            fd_para.row(i) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i - 1)) * 
                            (mK(i) / mK(i - 1) + mK(i - 1) / mK(i) * mH.row(i) / mH.row(i - 1)));
    }

    for (Eigen::Index i = mDipoleLayer; i < mNumLayers - 1; ++i) {
        c.row(i + 1) = 0.5 * Eigen::exp(I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                    (c.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (1 + mH.row(i) / mH.row(i + 1)) + 
                    (boolValue(i) + cd.row(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (1 - mH.row(i) / mH.row(i + 1)));

        cd.row(i + 1) = 0.5 * Eigen::exp(-I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                    (c.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (1 - mH.row(i) / mH.row(i + 1)) + 
                    (boolValue(i) + cd.row(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (1 + mH.row(i) / mH.row(i + 1)));

        f_perp.row(i + 1) = 0.5 * Eigen::exp(I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                            (f_perp.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) + mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)) + 
                            (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) - mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)));

        fd_perp.row(i + 1) = 0.5 * Eigen::exp(-I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                            (f_perp.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) - mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)) + 
                            (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) + mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)));

        f_para.row(i + 1) = 0.5 * Eigen::exp(I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                    (f_para.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (mK(i) / mK(i + 1) + mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)) + 
                    (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                    (mK(i) / mK(i + 1) - mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)));

        fd_para.row(i + 1) = 0.5 * Eigen::exp(-I * mH.row(i + 1) * (mZ0.cast<CMPLX>())(i)) * 
                            (f_para.row(i) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) - mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)) + 
                            (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)) * 
                            (mK(i) / mK(i + 1) + mK(i + 1) / mK(i) * mH.row(i) / mH.row(i + 1)));
    }

    mPowerPerpU.resize(mNumLayers - 1, mU.size() - 1);
    mPowerParaU.resize(mNumLayers - 1, mU.size() - 1);
    Eigen::ArrayXXcd powerParaTemp(mNumLayers - 1, mU.size() - 1);

    for (Eigen::Index i=0; i<mNumLayers-1; ++i) {
        mPowerPerpU.row(i) = (-3.0 * q * mdU / 4.0) * ((Eigen::pow(mU.segment(0, mU.size() - 1), 3)) / Eigen::abs(1 - Eigen::pow(mU.segment(0, mU.size() - 1), 2))) * 
                            (Eigen::sqrt(mEpsilon(i) / mEpsilon(mDipoleLayer) - Eigen::pow(mU.segment(0, mU.size() - 1), 2))) * 
                            (std::conj(std::sqrt(mEpsilon(i))) / std::sqrt(mEpsilon(i))); 
        mPowerPerpU.row(i) *= (f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) - 
                            ((fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)));
        mPowerPerpU.row(i) *= (Eigen::conj((f_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) +
                            ((fd_perp(i, Eigen::seqN(0, mPowerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)))));         
                            
        mPowerParaU.row(i) = (-3.0 * q * mdU / 8.0) * (mU.segment(0, mU.size() - 1) * Eigen::conj(Eigen::sqrt(mEpsilon(i) / mEpsilon(mDipoleLayer) - Eigen::pow(mU.segment(0, mU.size() - 1), 2)))) /
                            (Eigen::abs(1 - Eigen::pow(mU.segment(0, mU.size() - 1), 2))); 
        mPowerParaU.row(i) *= (c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) +
                            ((cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)));
        mPowerParaU.row(i) *= (Eigen::conj((c(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) -
                            ((cd(i, Eigen::seqN(0, mPowerParaU.cols())) + boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)))));
                            
        powerParaTemp.row(i) = (-3.0 * q * mdU / 8.0) * (mU.segment(0, mU.size() - 1) * Eigen::sqrt(mEpsilon(i)/mEpsilon(mDipoleLayer) - Eigen::pow(mU.segment(0, mU.size() - 1), 2))) *
                            (std::conj(std::sqrt(mEpsilon(i))) / std::sqrt(mEpsilon(i)));
        powerParaTemp.row(i) *= (f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) -
                                ((fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i)));
        powerParaTemp.row(i) *= (Eigen::conj((f_para(i, Eigen::seqN(0, mPowerParaU.cols())) * Eigen::exp(-I * mH.row(i) * (mZ0.cast<CMPLX>())(i))) +
                                ((fd_para(i, Eigen::seqN(0, mPowerParaU.cols())) - boolValue(i)) * Eigen::exp(I * mH.row(i) * (mZ0.cast<CMPLX>())(i))))); 
    }
    mPowerParaU += powerParaTemp;
}

void Solver::calculateEmissionSubstrate(Eigen::ArrayXd& thetaGlass, Eigen::ArrayXd& powerGlass) {
    double uCriticalGlass = std::real(std::sqrt(mEpsilon(mNumLayers - 1) / mEpsilon(mDipoleLayer)));
    auto uGlassIt = std::find_if(mU.begin(), mU.end(), [uCriticalGlass](auto a){ return a > uCriticalGlass; });
    auto uGlassIndex = uGlassIt - mU.begin();

    thetaGlass = Eigen::real(Eigen::acos(Eigen::sqrt(1 - mEpsilon(mDipoleLayer) / mEpsilon(mNumLayers - 1) * Eigen::pow(mU(Eigen::seq(1, uGlassIndex)), 2))));
    powerGlass = ((Eigen::real(mPowerPerpU(mNumLayers-2, Eigen::seq(1, uGlassIndex)))) * std::sqrt(std::real(mEpsilon(mNumLayers-1) / mEpsilon(mDipoleLayer))));
    powerGlass /= Eigen::tan(thetaGlass);
}



