#pragma once

#include <filesystem>

class MFile
{
private:
  std::filesystem::path mPath{};
  std::map<double, std::complex<double>> mSimulationData{};
  std::map<double, double> mFittingData{}; 
  const char mDelimiter = '\t'; 

  int load(bool mode=0);

public:
  explicit MFile(const std::filesystem::path& path);
  MFile(const std::filesystem::path& path, const char delimiter);

  const std::map<double, std::complex<double>>& getData() const;
};