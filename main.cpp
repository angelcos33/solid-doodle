#include <filesystem>
#include <iostream>

#include <linalg.hpp>
#include <material.hpp>

int main(int argc, char* argv[])
{
  // Test input data
  std::filesystem::path dataFile("./mat/alq3_literature.dat");

  Material alq3 = Material(dataFile, ',');
  std::vector<double> alq3Data = alq3.getRefIndex();

  // Test slicing. alq3Data contains (wvl, n, k) as | wvl | n | k | wvl | n | k | ...
  // Here we print the first 10 values of wavelength and n
  std::vector<double> wvl = slice(alq3Data, 0, 30, 3);
  print(wvl);
  std::vector<double> alq3N = slice(alq3Data, 1, 30, 3);
  print(alq3N);
}