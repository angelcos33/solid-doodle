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
  explicit Material(double realRefIndex, double imagRefIndex = 0.0);
  explicit Material(double wavelength, double realRefIndex, double imagRefIndex = 0.0);
  explicit Material(const std::filesystem::path& path, const char delimiter = '\t');

  std::complex<double> getRefIndex(double wavelength) const;
  std::complex<double> getEpsilon(double wavelength) const;
};