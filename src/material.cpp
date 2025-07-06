#include "material.hpp"
#include "files.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

Material::Material(double realRefIndex, double imagRefIndex)
{
  mRefIndices.insert(std::pair<double, std::complex<double>>{0.0, std::complex<double>(realRefIndex, imagRefIndex)});
}

Material::Material(double wavelength, double realRefIndex, double imagRefIndex)
{
  mRefIndices.insert(std::pair{wavelength, std::complex<double>(realRefIndex, imagRefIndex)});
}

Material::Material(const std::filesystem::path& path, const char delimiter)
{
  MFile mFile = MFile(path, delimiter);
  mRefIndices = mFile.getData();
}

std::complex<double> Material::getRefIndex(double wavelength) {

    std::complex<double> res;
    try {
        res = mRefIndices.at(wavelength);
        return res;
    }
    catch (const std::out_of_range& oor) {
        auto const uBound = mRefIndices.upper_bound(wavelength);
        if (uBound != mRefIndices.end()) {
            auto const lBound = std::prev(uBound);
            res.real((lBound->second).real() + (((uBound->second).real() - (lBound->second).real()) * (wavelength - lBound->first) / (uBound->first - lBound->first))); // Interpolation
            res.imag((lBound->second).imag() + (((uBound->second).imag() - (lBound->second).imag()) * (wavelength - lBound->first) / (uBound->first - lBound->first))); // Interpolation
            return res;
        }
        else throw std::runtime_error("The wavelength provided is outside the range of the material data");
    }
}

std::complex<double> Material::getEpsilon(double wavelength)
{
  double prevKey = 0.0;
  std::complex<double> prevVal = 0.0;
  std::complex<double> epsilon = 0.0;
  ;

  for (std::map<double, std::complex<double>>::iterator iter = mRefIndices.begin(); iter != mRefIndices.end(); iter++) {
    double key = iter->first;
    std::complex<double> val = iter->second;

    double ratio = (key - wavelength) / (key - prevKey);
    if (key > wavelength && prevKey != 0.0) {
      epsilon = ratio * prevVal + (1 - ratio) * val;
      epsilon = std::pow(epsilon, 2);
      return epsilon;
    }
    else if (key < wavelength && std::next(iter, 1) != mRefIndices.end()) {
      epsilon = std::pow(val, 2);
      std::cerr << "WARNING: Wavelength selected is below the range available, epsilon might be wrong \n";
      return epsilon;
    }
    else if (key > wavelength) {
      epsilon = std::pow(val, 2);
      std::cerr << "WARNING: Wavelength selected is below the range available, epsilon might be wrong \n";
      return epsilon;
    }
    else if (key == wavelength) { // exact value
      epsilon = std::pow(val, 2);
      return epsilon;
    }
    prevKey = key;
    prevVal = val;
  }
}