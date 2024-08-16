/*
Generic library for vector utilities. It is based on the std::vector container.
*/
#pragma once

#include<iostream>
#include<vector>
#include<cassert>

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
    for (size_t i=0; i<vec.size(); ++i) {
        std::cout << vec[i] << "\n";
    }
}