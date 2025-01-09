#include<iostream>

#include<./linalg.hpp>

int main(int argc, char* argv[]) {
    std::vector<double> v;
    linspace(v, 1.1, 2.0, 10.0);
    std::vector<double> u (10, 1.0);
    print(v);
    std::cout << "\n";
    // Test vector multiplication
    std::vector<double> r = v * u;
    print(r);
    std::cout << "\n";
    print(u);
}