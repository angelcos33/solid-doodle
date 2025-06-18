#pragma once

#include <complex>
#include <filesystem>
#include <map>
#include <vector>

class Material
{
private:
  std::map<double, std::complex<double>> mRefIndices;

public:
  Material(double realRefIndex, double imagRefIndex = 0.0);
  Material(double wavelength, double realRefIndex, double imagRefIndex = 0.0);
  Material(const std::filesystem::path& path, const char delimiter = '\t');

  std::complex<double> getRefIndex(double wavelength);
  std::complex<double> getEpsilon(double wavelength);
};