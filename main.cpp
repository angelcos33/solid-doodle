#include<iostream>
#include<filesystem>

#include<linalg.hpp>
#include<files.hpp>
#include<material.hpp>

int main(int argc, char* argv[]) {
    // Test input data
    std::filesystem::path dataFile ("./mat/alq3_literature.dat");


    MFile mfile (dataFile, ',');
    std::vector<double> data = mfile.getData();
    
    Material alq3 = Material(data);
}