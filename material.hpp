#pragma once

#include<vector>

class Material {
    private:
        std::vector<double> mRefIndex{};

    public:
        explicit Material(double realRefIndex);
        Material(double realRefIndex, double imagRefIndex);
        Material(std::vector<double> refIndexData);

        std::vector<double> getRefIndex() const;
};