/*
Generic library for vector utilities. It is based on the std::vector container.
*/
#pragma once

#include<iostream>
#include<vector>
#include<cassert>

template<typename T>
void linspace(std::vector<T>& out, double start, double stop, size_t N) {
    out.reserve(N);

    assert(stop > start);
    double step = (stop - start) / (N - 1);

    for (size_t i=0; i<N; ++i) {
        out.push_back(start + i*step);
    }
}

template<typename T>
void print(std::vector<T> vec) {
    for (size_t i=0; i<vec.size(); ++i) {
        std::cout << vec[i] << "\n";
    }
}