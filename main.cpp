#include<iostream>

#include<./linalg.hpp>

int main(int argc, char* argv[]) {
    std::vector<double> v; 
    linspace(v, 1.0f, 10.0f, 10);
    print(v);
}