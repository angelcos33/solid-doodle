#pragma once

#include<vector>
#include<filesystem>

class Material {
    private:
        std::vector<double> mRefIndex{};

    public:
        explicit Material(double realRefIndex);
        Material(double realRefIndex, double imagRefIndex);
        Material(const std::filesystem::path& path, const char delimiter='\t');

        std::vector<double> getRefIndex() const;
};
