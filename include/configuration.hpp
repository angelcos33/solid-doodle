#pragma once

#include <vector>

#include <simulation.hpp>
#include <jsonsimplecpp/parser.hpp>

struct Input {
    std::vector<Material> materials;
    std::vector<double> thicknesses;
    double emitterPosition;
    Eigen::Index emitterIndex;
    double wavelength;
    GaussianSpectrum spectrum;
    DipoleDistribution dipoleDist;
};

class ConfigurationManager {
   Json::JsonParser jsonParser; 

   public:
    ConfigurationManager(const std::string& filename);

    void configure();

};