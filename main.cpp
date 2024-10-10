#include<iostream>
#include<filesystem>

#include<linalg.hpp>
#include<files.hpp>

int main(int argc, char* argv[]) {
    // Test input data
    std::filesystem::path dataFile ("./dummy_nk.nk");


    MFile mfile (dataFile);
    std::vector<double> data = mfile.getData();
    print(data);
}