#include<algorithm>
#include<material.hpp>
#include<files.hpp>

Material::Material(double realRefIndex) {
    mRefIndex.push_back(realRefIndex);
}

Material::Material(double realRefIndex, double imagRefIndex) {
    mRefIndex.push_back(realRefIndex);
    mRefIndex.push_back(imagRefIndex);
}

Material::Material(const std::filesystem::path& path, const char delimiter) {
    MFile mFile = MFile(path, delimiter);
    mRefIndex = mFile.getData();
}

std::vector<double> Material::getRefIndex() const{
    return mRefIndex;
}
