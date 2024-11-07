#include<iostream>
#include<filesystem>

#include<linalg.hpp>
#include<material.hpp>

int main(int argc, char* argv[]) {
    // Test input data
    std::filesystem::path dataFile ("./mat/alq3_literature.dat");

    Material alq3 = Material(dataFile, ',');
    print(alq3.getRefIndex());
}