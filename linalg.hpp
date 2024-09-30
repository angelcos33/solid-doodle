/*
Generic library for vector utilities. It is based on the std::vector container.
*/
#pragma once

#include<iostream>
#include<vector>
#include<cassert>

// Vector arithmetic
// Scalar Multiplication
template<typename T, typename U>
std::vector<T> operator*(const std::vector<T>& vec, const U scalar) {
    std::vector<T> res;
    res.reserve(vec.size());

    for (auto& elem: vec) {
        res.push_back(scalar * elem);
    }
    return res;
}

template<typename T, typename U>
std::vector<T> operator*(const U scalar, const std::vector<T>& vec) {
    std::vector<T> res;
    res.reserve(vec.size());

    for (auto& elem: vec) {
        res.push_back(scalar * elem);
    }
    return res;
}

// Vector addition
template<typename T>
std::vector<T> operator+(const std::vector<T>& lvec, const std::vector<T>& rvec) {
    assert(lvec.size() == rvec.size());
    std::vector<T> res;
    res.reserve(lvec.size());

    for (size_t i=0; i<lvec.size(); ++i) {
        res.push_back(lvec[i] + rvec[i]);
    }
    return res;
}

// Vector element-wise multiplication
template<typename T>
std::vector<T> operator*(const std::vector<T>& lvec, const std::vector<T>& rvec) {
    assert(lvec.size() == rvec.size());
    std::vector<T> res;
    res.reserve(lvec.size());

    for (size_t i=0; i<lvec.size(); ++i) {
        res.push_back(lvec[i] * rvec[i]);
    }
    return res;
}

template<typename T>
void arange(std::vector<T>& out, T start, T stop, T step) {

    assert(stop > start);
    size_t N = ((stop - start) / step) + 1;
    out.reserve(N);

    for (size_t i=0; i<N; ++i) {
        out.push_back(start + i*step);
    }
}

template<typename T>
void linspace(std::vector<T>& out, T start, T stop, size_t N) {
    T step = (stop - start) / (N - 1);
    arange(out, start, stop, step);
}

template<typename T>
void print(const std::vector<T>& vec) {
    for (auto& elem: vec) {
        std::cout << elem << "\n";
    }
}