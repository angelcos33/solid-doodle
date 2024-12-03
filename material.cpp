#include<algorithm>
#include<iostream>
#include<complex>
#include<cmath>
#include"material.hpp"
#include"files.hpp"

Material::Material(double realRefIndex, double imagRefIndex = 0.0) {
    mRefIndices.insert(std::pair<double, std::complex<double>>{0.0, std::complex<double>(realRefIndex, imagRefIndex)});
}

Material::Material(double wavelength, double realRefIndex, double imagRefIndex=0.0) {
    mRefIndices.insert(std::pair{wavelength, std::complex<double>(realRefIndex, imagRefIndex)});
}

Material::Material(const std::filesystem::path& path, const char delimiter) {
    MFile mFile = MFile(path, delimiter);
    mRefIndices = mFile.getData();
}

std::complex<double> Material::getRefIndex(double wavelength) {
    double prevKey = 0.0;
    std::complex<double> prevVal = 0.0;
    std::complex<double> refIndex = 0.0;

    for(std::map<double, std::complex<double>>::iterator iter = mRefIndices.begin(); iter!=mRefIndices.end(); iter++) {
        double key = iter->first;
        std::complex<double> val = iter->second;

        double ratio = (key-wavelength)/(key-prevKey);
        if(key > wavelength && prevKey != 0.0) {
            refIndex = ratio*prevVal + (1-ratio)*val;
            return refIndex;
        }
        else if(key < wavelength && std::next(iter, 1) != mRefIndices.end()) {
            refIndex = val;
            std::cerr << "WARNING: Wavelength selected is above the range available, n might be wrong \n"; 
            return refIndex;
        }
        else if(key > wavelength) {
            refIndex = val;
            std::cerr << "WARNING: Wavelength selected is below the range available, n might be wrong \n";
            return refIndex;
        }
        else if (key == wavelength) { //exact value
            refIndex = val;
            return refIndex;
        }
        prevKey = key;
        prevVal = val;
    }
}

std::complex<double> Material::getEpsilon(double wavelength) {
    double prevKey = 0.0;
    std::complex<double> prevVal = 0.0;
    std::complex<double> epsilon = 0.0;;

    for(std::map<double, std::complex<double>>::iterator iter = mRefIndices.begin(); iter!=mRefIndices.end(); iter++) {
        double key = iter->first;
        std::complex<double> val = iter->second;

        double ratio = (key-wavelength)/(key-prevKey);
        if(key > wavelength && prevKey != 0.0) {
            epsilon = ratio*prevVal + (1-ratio)*val;
            epsilon = std::pow(epsilon, 2);
            return epsilon;
        }
        else if(key < wavelength && std::next(iter, 1) != mRefIndices.end()) {
            epsilon = std::pow(val, 2);
            std::cerr << "WARNING: Wavelength selected is below the range available, epsilon might be wrong \n";
            return epsilon;
        }
        else if(key > wavelength) {
            epsilon = std::pow(val, 2);
            std::cerr << "WARNING: Wavelength selected is below the range available, epsilon might be wrong \n";
            return epsilon;
        } 
        else if (key == wavelength) { //exact value
            epsilon = std::pow(val, 2);
            return epsilon;
        }
        prevKey = key;
        prevVal = val;
    }
} 