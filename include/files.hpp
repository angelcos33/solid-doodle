/*! \file files.hpp
    \brief a header file which contains a small class to interface with files.

    The files.hpp header was designed to contained all the needed file input/output 
    functionality. This entails primarily the ability to read experimental data and 
    material information.
*/

#pragma once

#include <filesystem>

/*! \class MFile
    \brief a simple file parsing and information extraction class.
    
    The MFile class provides convinient file parsing and data extraction
    finctionality. Thus, it acts as the interface between the 
    raw material/experimental data and the solver classes used for
    simulation and fitting.
*/ 
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