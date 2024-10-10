#include<filesystem>
#include<iostream>
#include<cassert>
#include<fstream>
#include<string>
#include<sstream>
#include<cctype>
#include<algorithm>

#include<./files.hpp>

MFile::MFile(const std::filesystem::path& path): mPath{path} {
    this->load();
}

std::vector<double> MFile::getData() const{
    return mData;
}

int MFile::load() {
    assert(std::filesystem::exists(mPath));

    std::ifstream iFile(mPath);

    if (!iFile.is_open()) {
        std::cerr << "Failed to open file " << mPath << "\n";
        return 1;
    }

    std::string line;
    std::vector<size_t> cols;
    bool startOfDataFound = false;

    while(std::getline(iFile, line)) {
        // Skip header lines, i.e. if they start with non digit values
        if (!(startOfDataFound)) {
            if (!(line.empty()) && std::isdigit(line[0])) {
                startOfDataFound = true;
            }
            else continue;
        }
        std::istringstream iss(line);
        std::string token;
        size_t col = 0;
        while(std::getline(iss, token, mDelimiter)) {
                try
                {   
                    col++;
                    mData.push_back(std::stod(token));
                }
                catch(const std::invalid_argument& e)
                {
                    std::cerr << "Conversion error: " << e.what() << "for token: " << token << '\n';
                }
                catch(const std::out_of_range& e) {
                    std::cerr << "Out of range error: " << e.what() << "for token: " << token << "\n";
                }
        }
        cols.push_back(col);
    }
    // Make sure that all the columns have same number of values!
    assert(std::all_of(cols.begin(), cols.end(),
                        [first = cols.front()](size_t value) {return value == first;}));
    // Also make sure there are only 2 columns, 1 for n and 1 for k!
    assert(std::all_of(cols.begin(), cols.end(),
                        [](size_t value) {return value == 2;}));                        
    iFile.close();

}