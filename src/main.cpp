#include <iostream>
#include <linalg.hpp>
#include <matplotlib.hpp>
#include <algorithm>
#include <complex>
#include <cmath>
#include <iterator>
#include <numeric>
#include <functional>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

using CMPLX = std::complex<double>;

int main()
{
  /*
  std::vector<Material> materials; 
  std::vector<double> thicknesses;
  // Air
  materials.push_back(Material(1.0));
  thicknesses.push_back(1.0);
  // Film
  materials.push_back(Material(2));
  thicknesses.push_back(50e-9);
  // Glass
  materials.push_back(Material(1.5));
  thicknesses.push_back(1.0);

  // Create MatStack
  MatStack stack = MatStack(1,
                            500e-9,
                            thicknesses,
                            materials);
  size_t N = 1000000;
  
  ArrayXd a = ArrayXd::Random(N);
  ArrayXd res(N);

  auto start = std::chrono::steady_clock::now();
  std::partial_sum(a.begin(), a.end(), res.begin(), std::plus<double>());
  auto stop = std::chrono::steady_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
  fmt::print("Runtime of cumsum: {:4}\n", elapsed);
  */
  Eigen::Index N = 10; // Num of layers
  double q = 1.0; // PLQY
  Eigen::Index dipole_layer = 4;
  double dipole_position = 0.0;
  Eigen::ArrayXd thickness(N - 2), z0(N - 1);
  Eigen::ArrayXcd epsilon(N);
  thickness << 200e-10, 1000e-10, 500e-10, 200e-10, 500e-10, 300e-10, 1600e-10, 5000e-10;
  epsilon << CMPLX(1.0, 0.0),
             CMPLX(-14.924227573900003, 0.785503814000000),
             CMPLX(-24.819095227053140, 3.137674975845411),
             CMPLX(2.985292929845569, 0.000041123270350),
             CMPLX(2.961496810000000, 0.0),
             CMPLX(2.985292929845569, 0.000041123270350),
             CMPLX(2.304131952206292, 0.033315508532326),
             CMPLX(3.294650000000000, 0.036300000000000),
             CMPLX(2.25, 0.0),
             CMPLX(2.25, 0.0);
  // Cumulative sum
  z0(0) = 0.0;
  std::partial_sum(thickness.begin(), thickness.end(), std::next(z0.begin()), std::plus<double>());
  z0 -= (z0(dipole_layer - 1) + dipole_position);
  // Discretization
  CMPLX I(0.0, 1.0);
  double x_res = 5e-4;
  Eigen::ArrayXd x_real = arange<Eigen::ArrayXd>(-M_PI_2, -x_res, x_res);
  Eigen::ArrayXcd x_imag = I * arange<Eigen::ArrayXd>(x_res, 1.8, x_res);
  Eigen::ArrayXcd x(x_real.rows() + x_imag.rows());
  x.head(x_real.size()) = x_real.cast<CMPLX>();
  x.segment(x_real.size(), x_imag.size()) = x_imag;
  Eigen::ArrayXd u = x.cos().real();
  // Differences
  Eigen::ArrayXcd du = u.segment(1, u.size() - 1) - u.segment(0, u.size() - 1);

  double wavelength = 535; // nm
  Eigen::ArrayXcd k = 2 * M_PI / wavelength / 1e-9 * epsilon.sqrt();

  // h: Out of plane normalized wavevector k_z
  Eigen::ArrayXXcd h(epsilon.size(), x.size());
  h = k(dipole_layer) * (((epsilon.replicate(1, x.size()))/epsilon(dipole_layer)).rowwise() - x.cos().pow(2).transpose()).sqrt();

  // R_para/perp: Fresnel coefficients
  Eigen::ArrayXXcd R_perp(epsilon.size() - 1, x.size()), R_para(epsilon.size() - 1, x.size());
  R_perp = (h.block(1, 0, h.rows() - 1, h.cols()) - h.block(0, 0, h.rows() - 1, h.cols())) / 
           (h.block(1, 0, h.rows() - 1, h.cols()) + h.block(0, 0, h.rows() - 1, h.cols()));
  (R_perp.bottomRows(R_perp.rows() - dipole_layer)) *= -1.0; 

  R_para = ((h.block(0, 0, h.rows() - 1, h.cols())).colwise() * epsilon.segment(1, epsilon.size() - 1) - 
            (h.block(1, 0, h.rows() - 1, h.cols())).colwise() * epsilon.segment(0, epsilon.size() - 1));
  R_para /= ((h.block(0, 0, h.rows() - 1, h.cols())).colwise() * epsilon.segment(1, epsilon.size() - 1) + 
            (h.block(1, 0, h.rows() - 1, h.cols())).colwise() * epsilon.segment(0, epsilon.size() - 1)); 
  (R_para.bottomRows(R_para.rows() - dipole_layer)) *= -1.0; 

  // CB and FB coefficients
  Eigen::ArrayXXcd CB = Eigen::ArrayXXcd::Zero(dipole_layer + 1, x.size());
  Eigen::ArrayXXcd FB = Eigen::ArrayXXcd::Zero(dipole_layer + 1, x.size());
  for (Eigen::Index i = 1; i < dipole_layer + 1; ++i) {
    CB.row(i) = Eigen::exp(-2.0 * I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * (R_perp.row(i-1) + (CB.row(i-1) * Eigen::exp(2.0 * I * h.row(i-1) * (z0.cast<CMPLX>())(i - 1)))) / 
                (1 + R_perp.row(i - 1) * (CB.row(i-1) * Eigen::exp(2.0 * I * h.row(i-1) * (z0.cast<CMPLX>())(i - 1))));

    FB.row(i) = Eigen::exp(-2.0 * I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * (-R_para.row(i-1) + (FB.row(i-1) * Eigen::exp(2.0 * I * h.row(i-1) * (z0.cast<CMPLX>())(i - 1)))) / 
                (1 - R_para.row(i - 1) * (FB.row(i-1) * Eigen::exp(2.0 * I * h.row(i-1) * (z0.cast<CMPLX>())(i - 1))));
  }

  // CT and FT coefficients
  Eigen::ArrayXXcd CT = Eigen::ArrayXXcd::Zero(N - dipole_layer, x.size());
  Eigen::ArrayXXcd FT = Eigen::ArrayXXcd::Zero(N - dipole_layer, x.size());
  for (Eigen::Index i = N - dipole_layer - 2; i >= 0; --i) {
    Eigen::Index indexFromTop = i + dipole_layer;
    CT.row(i) = Eigen::exp(2.0 * I * h.row(indexFromTop) * (z0.cast<CMPLX>())(indexFromTop)) * (R_perp.row(indexFromTop) + (CT.row(i+1) * Eigen::exp(-2.0 * I * h.row(indexFromTop+1) * (z0.cast<CMPLX>())(indexFromTop)))) / 
                (1 + R_perp.row(indexFromTop) * (CT.row(i+1) * Eigen::exp(-2.0 * I * h.row(indexFromTop + 1) * (z0.cast<CMPLX>())(indexFromTop))));

    FT.row(i) = Eigen::exp(2.0 * I * h.row(indexFromTop) * (z0.cast<CMPLX>())(indexFromTop)) * (-R_para.row(indexFromTop) + (FT.row(i+1) * Eigen::exp(-2.0 * I * h.row(indexFromTop+1) * (z0.cast<CMPLX>())(indexFromTop)))) / 
                (1 - R_para.row(indexFromTop) * (FT.row(i+1) * Eigen::exp(-2.0 * I * h.row(indexFromTop + 1) * (z0.cast<CMPLX>())(indexFromTop))));
  }

  // c, cd, f, fd coefficients
  Eigen::ArrayXXcd c = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());
  Eigen::ArrayXXcd cd = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());
  Eigen::ArrayXXcd f_perp = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());
  Eigen::ArrayXXcd fd_perp = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());
  Eigen::ArrayXXcd f_para = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());
  Eigen::ArrayXXcd fd_para = Eigen::ArrayXXcd::Zero(epsilon.size(), x.size());

  c.row(dipole_layer) = (CB.row(dipole_layer) + 1) * CT.row(0) / (1 - CB.row(dipole_layer) * CT.row(0));
  cd.row(dipole_layer) = (CT.row(0) + 1) * CB.row(dipole_layer) / (1 - CB.row(dipole_layer) * CT.row(0));

  f_perp.row(dipole_layer) = (FB.row(dipole_layer) + 1) * FT.row(0) / (1 - FB.row(dipole_layer) * FT.row(0));
  fd_perp.row(dipole_layer) = (FT.row(0) + 1) * FB.row(dipole_layer) / (1 - FB.row(dipole_layer) * FT.row(0));

  f_para.row(dipole_layer) = (FB.row(dipole_layer) - 1) * FT.row(0) / (1 - FB.row(dipole_layer) * FT.row(0));
  fd_para.row(dipole_layer) = (1 - FT.row(0)) * FB.row(dipole_layer) / (1 - FB.row(dipole_layer) * FT.row(0));

  Eigen::ArrayXd boolValue = Eigen::ArrayXd::Zero(epsilon.size());
  boolValue(dipole_layer) = 1.0;

  for (Eigen::Index i = dipole_layer; i>=1; --i) {
    c.row(i - 1) = 0.5 * Eigen::exp(I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                   ((boolValue(i) + c.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (1 + h.row(i) / h.row(i - 1)) + 
                   cd.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (1 - h.row(i) / h.row(i - 1)));

    cd.row(i - 1) = 0.5 * Eigen::exp(-I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                   ((boolValue(i) + c.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (1 - h.row(i) / h.row(i - 1)) + 
                   cd.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (1 + h.row(i) / h.row(i - 1)));

    f_perp.row(i - 1) = 0.5 * Eigen::exp(I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                        ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                        (k(i) / k(i - 1) + k(i - 1) / k(i) * h.row(i) / h.row(i - 1)) + 
                        fd_perp.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                        (k(i) / k(i - 1) - k(i - 1) / k(i) * h.row(i) / h.row(i - 1)));

    fd_perp.row(i - 1) = 0.5 * Eigen::exp(-I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                         ((boolValue(i) + f_perp.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                         (k(i) / k(i - 1) - k(i - 1) / k(i) * h.row(i) / h.row(i - 1)) + 
                         fd_perp.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                         (k(i) / k(i - 1) + k(i - 1) / k(i) * h.row(i) / h.row(i - 1)));

    f_para.row(i - 1) = 0.5 * Eigen::exp(I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                   ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (k(i) / k(i - 1) + k(i - 1) / k(i) * h.row(i) / h.row(i - 1)) + 
                   fd_para.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                   (k(i) / k(i - 1) - k(i - 1) / k(i) * h.row(i) / h.row(i - 1)));

    fd_para.row(i - 1) = 0.5 * Eigen::exp(-I * h.row(i - 1) * (z0.cast<CMPLX>())(i - 1)) * 
                         ((boolValue(i) + f_para.row(i)) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                         (k(i) / k(i - 1) - k(i - 1) / k(i) * h.row(i) / h.row(i - 1)) + 
                         fd_para.row(i) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i - 1)) * 
                         (k(i) / k(i - 1) + k(i - 1) / k(i) * h.row(i) / h.row(i - 1)));
  }

  for (Eigen::Index i = dipole_layer; i < epsilon.size()-1; ++i) {
    c.row(i + 1) = 0.5 * Eigen::exp(I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                   (c.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (1 + h.row(i) / h.row(i + 1)) + 
                   (boolValue(i) + cd.row(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (1 - h.row(i) / h.row(i + 1)));

    cd.row(i + 1) = 0.5 * Eigen::exp(-I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                   (c.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (1 - h.row(i) / h.row(i + 1)) + 
                   (boolValue(i) + cd.row(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (1 + h.row(i) / h.row(i + 1)));

    f_perp.row(i + 1) = 0.5 * Eigen::exp(I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                        (f_perp.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                        (k(i) / k(i + 1) + k(i + 1) / k(i) * h.row(i) / h.row(i + 1)) + 
                        (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                        (k(i) / k(i + 1) - k(i + 1) / k(i) * h.row(i) / h.row(i + 1)));

    fd_perp.row(i + 1) = 0.5 * Eigen::exp(-I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                         (f_perp.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                         (k(i) / k(i + 1) - k(i + 1) / k(i) * h.row(i) / h.row(i + 1)) + 
                         (boolValue(i) + fd_perp.row(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                         (k(i) / k(i + 1) + k(i + 1) / k(i) * h.row(i) / h.row(i + 1)));

    f_para.row(i + 1) = 0.5 * Eigen::exp(I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                   (f_para.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (k(i) / k(i + 1) + k(i + 1) / k(i) * h.row(i) / h.row(i + 1)) + 
                   (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                   (k(i) / k(i + 1) - k(i + 1) / k(i) * h.row(i) / h.row(i + 1)));

    fd_para.row(i + 1) = 0.5 * Eigen::exp(-I * h.row(i + 1) * (z0.cast<CMPLX>())(i)) * 
                         (f_para.row(i) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                         (k(i) / k(i + 1) - k(i + 1) / k(i) * h.row(i) / h.row(i + 1)) + 
                         (fd_para.row(i) - boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)) * 
                         (k(i) / k(i + 1) + k(i + 1) / k(i) * h.row(i) / h.row(i + 1)));
  }

  // Poynting vector at interfaces
  Eigen::ArrayXXcd powerPerpU(epsilon.size() - 1, x.size() - 1);
  Eigen::ArrayXXcd powerParaU(epsilon.size() - 1, x.size() - 1);
  Eigen::ArrayXXcd powerParaTemp(epsilon.size() - 1, x.size() - 1);

  
  for (Eigen::Index i=0; i<epsilon.size()-1; ++i) {
    powerPerpU.row(i) = (-3.0 * q * du / 4.0) * ((Eigen::pow(u.segment(0, u.size() - 1), 3)) / Eigen::abs(1 - Eigen::pow(u.segment(0, u.size() - 1), 2))) * 
                        (Eigen::sqrt(epsilon(i) / epsilon(dipole_layer) - Eigen::pow(u.segment(0, u.size() - 1), 2))) * 
                        (std::conj(std::sqrt(epsilon(i))) / std::sqrt(epsilon(i))); 
    powerPerpU.row(i) *= (f_perp(i, Eigen::seqN(0, powerPerpU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) - 
                         ((fd_perp(i, Eigen::seqN(0, powerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)));
    powerPerpU.row(i) *= (Eigen::conj((f_perp(i, Eigen::seqN(0, powerPerpU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) +
                         ((fd_perp(i, Eigen::seqN(0, powerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)))));         
                         
    powerParaU.row(i) = (-3.0 * q * du / 8.0) * (u.segment(0, u.size() - 1) * Eigen::conj(Eigen::sqrt(epsilon(i) / epsilon(dipole_layer) - Eigen::pow(u.segment(0, u.size() - 1), 2)))) /
                        (Eigen::abs(1 - Eigen::pow(u.segment(0, u.size() - 1), 2))); 
    powerParaU.row(i) *= (c(i, Eigen::seqN(0, powerParaU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) +
                         ((cd(i, Eigen::seqN(0, powerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)));
    powerParaU.row(i) *= (Eigen::conj((c(i, Eigen::seqN(0, powerPerpU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) -
                         ((cd(i, Eigen::seqN(0, powerPerpU.cols())) + boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)))));
                        
    powerParaTemp.row(i) = (-3.0 * q * du / 8.0) * (u.segment(0, u.size() - 1) * Eigen::sqrt(epsilon(i)/epsilon(dipole_layer) - Eigen::pow(u.segment(0, u.size() - 1), 2))) *
                           (std::conj(std::sqrt(epsilon(i))) / std::sqrt(epsilon(i)));
    powerParaTemp.row(i) *= (f_para(i, Eigen::seqN(0, powerParaU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) -
                            ((fd_para(i, Eigen::seqN(0, powerParaU.cols())) - boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i)));
    powerParaTemp.row(i) *= (Eigen::conj((f_para(i, Eigen::seqN(0, powerParaU.cols())) * Eigen::exp(-I * h.row(i) * (z0.cast<CMPLX>())(i))) +
                            ((fd_para(i, Eigen::seqN(0, powerPerpU.cols())) - boolValue(i)) * Eigen::exp(I * h.row(i) * (z0.cast<CMPLX>())(i))))); 
  }
  powerParaU += powerParaTemp;

  // Polar figure
  double uCriticalGlass = std::real(std::sqrt(epsilon(N - 1) / epsilon(dipole_layer)));
  auto uGlassIt = std::find_if(u.begin(), u.end(), [uCriticalGlass](auto a){ return a > uCriticalGlass; });
  auto uGlassIndex = uGlassIt - u.begin();
  //p_angle_glass=(power_perp_u(t,[2:u_glass],N-1).*sqrt(epsilon(N)/epsilon(dipole_layer))./tan(theta_glass));
  Eigen::ArrayXd thetaGlass = Eigen::real(Eigen::acos(Eigen::sqrt(1 - epsilon(dipole_layer) / epsilon(N - 1) * Eigen::pow(u(Eigen::seq(1, uGlassIndex)), 2))));
  Eigen::ArrayXd powerAngleGlass = ((Eigen::real(powerPerpU(N-2, Eigen::seq(1, uGlassIndex)))) * std::sqrt(std::real(epsilon(N-1) / epsilon(dipole_layer))));
  powerAngleGlass /= Eigen::tan(thetaGlass);
  std::cout << thetaGlass.head(10) << '\n';
  std::cout << powerAngleGlass.head(10) << '\n';

  figure();
  plot(thetaGlass, powerAngleGlass);
  show();
}