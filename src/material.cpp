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

std::complex<double> Material::getRefIndex(double wavelength) const
{

  std::complex<double> res;
  try {
    res = mRefIndices.at(wavelength);
    return res;
  } catch (const std::out_of_range& oor) {
    auto const uBound = mRefIndices.upper_bound(wavelength);
    if (uBound != mRefIndices.end()) {
      auto const lBound = std::prev(uBound);
      double ratio = (wavelength - lBound->first) / (uBound->first - lBound->first);
      res = ratio * uBound->second + (1 - ratio) * lBound->second; // Interpolation
      return res;
    }
    else throw std::runtime_error("The wavelength provided is outside the range of the material data");
  }
}

std::complex<double> Material::getEpsilon(double wavelength) const
{
  std::complex<double> res = getRefIndex(wavelength);
  return res * res;
}