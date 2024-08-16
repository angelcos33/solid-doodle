#include<iostream>

#include<./linalg.hpp>

int main(int argc, char* argv[]) {
    std::vector<double> v;
    linspace(v, 1.1, 2.0, 10.0);
    print(v);
}