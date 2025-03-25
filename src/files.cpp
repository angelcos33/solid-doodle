#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <files.hpp>
#include <massert.hpp>

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

std::vector<double> MFile::getData() const { return mData; }

int MFile::load()
{
  m_assert(std::filesystem::exists(mPath), "File doesn't exist!");

  std::ifstream iFile(mPath);

  if (!iFile.is_open()) {
    std::cerr << "Failed to open file " << mPath << "\n";
    return 1;
  }

  std::string line;
  std::vector<size_t> cols;
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
        col++;
        mData.push_back(std::stod(token));
      } catch (const std::invalid_argument& e) {
        std::cerr << "Conversion error: " << e.what() << "for token: " << token << '\n';
      } catch (const std::out_of_range& e) {
        std::cerr << "Out of range error: " << e.what() << "for token: " << token << "\n";
      }
    }
    cols.push_back(col);
  }
  // Make sure that all the columns have same number of values!
  assert(std::all_of(cols.begin(), cols.end(), [first = cols.front()](size_t value) { return value == first; }));
  // Also make sure there are only 3 columns, 1 for wvl, 2 for n and 3 for k!
  assert(std::all_of(cols.begin(), cols.end(), [](size_t value) { return value == MAX_COLS; }));
  iFile.close();
  return 0;
}