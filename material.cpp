#include<algorithm>
#include<./material.hpp>

Material::Material(double realRefIndex) {
    mRefIndex.push_back(realRefIndex);
}

Material::Material(double realRefIndex, double imagRefIndex) {
    mRefIndex.push_back(realRefIndex);
    mRefIndex.push_back(imagRefIndex);
}

Material::Material(std::vector<double> refIndexData) {
    mRefIndex = std::move(refIndexData);
}

std::vector<double> Material::getRefIndex() const{
    return mRefIndex;
}
