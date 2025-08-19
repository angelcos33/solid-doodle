#include <algorithm>
#include <array>
#include <cctype>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "files.hpp"
#include "massert.hpp"

#define MAX_COLS 3

MFile::MFile(const std::filesystem::path& path) :
  mPath{path}
{
  this->load();
}

MFile::MFile(const std::filesystem::path& path, const char delimiter) :
  mPath{path},
  mDelimiter{delimiter}
{
  this->load();
}

const std::map<double, std::complex<double>>& MFile::getData() const { return mData; }

int MFile::load()
{
  std::ifstream iFile(mPath);

  if (!iFile.is_open()) { std::perror(("Error while opening file " + mPath.string()).c_str()); }

  std::string line;
  std::vector<size_t> cols;
  std::array<double, MAX_COLS> data{};
  bool startOfDataFound = false;

  while (std::getline(iFile, line)) {
    // Skip header lines, i.e. if they start with non digit values
    if (!(startOfDataFound)) {
      if (!(line.empty()) && std::isdigit(line[0])) { startOfDataFound = true; }
      else continue;
    }
    std::istringstream iss(line);
    std::string token;
    size_t col = 0;
    while (std::getline(iss, token, mDelimiter)) {
      try {
        data[col] = (std::stod(token));
        col++;
      } catch (const std::invalid_argument& e) {
        std::cerr << "Conversion error: " << e.what() << "for token: " << token << '\n';
      } catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << "for token: " << token << "\n";
      }
    }
    cols.push_back(col);
    mData.insert(std::pair{data[MAX_COLS - 3], std::complex<double>{data[MAX_COLS - 2], data[MAX_COLS - 1]}});
  }
  if (iFile.bad()) std::perror(("Error while reading the file " + mPath.string()).c_str());
  // Make sure that all the columns have same number of values!
  m_assert(std::all_of(cols.begin(), cols.end(), [first = cols.front()](size_t value) { return value == first; }),
    "Material file cannot be parsed. Columns have different number of lines");
  // Also make sure there are only 3 columns, 1 for wvl, 2 for n and 3 for k!
  m_assert(std::all_of(cols.begin(), cols.end(), [](size_t value) { return value == MAX_COLS; }),
    "Material file cannot be parse. Number of columns must be 3");
  iFile.close();
  return 0;
}