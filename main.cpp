#include<iostream>

#include<./linalg.hpp>

int main(int argc, char* argv[]) {
    std::vector<double> v;
    linspace(v, 1.1, 2.0, 10.0);
    std::vector<double> u (10, 1.0);
    print(v);
    std::cout << "\n";
    // Test scalar multiplication
    std::vector<double> r = 2.0 * v;
    print(v);
    std::cout << "\n";
    print(r);
}