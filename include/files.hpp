#pragma once

#include <filesystem>

class MFile
{
private:
  bool mMode{};
  std::filesystem::path mPath{};
  std::map<double, std::complex<double>> mSimulationData{};
  std::map<double, double> mFittingData{}; 
  const char mDelimiter = '\t'; 

  int load();

public:
  explicit MFile(const std::filesystem::path& path, bool mode);
  MFile(const std::filesystem::path& path, const char delimiter, bool mode);

  void setNewPath(const std::filesystem::path& newPath);
  void setNewMode(bool newMode);

  const std::map<double, std::complex<double>>& getData() const;
};