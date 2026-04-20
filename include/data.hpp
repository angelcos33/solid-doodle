#pragma once

#include <string>

#include <forwardDecl.hpp>

class Data {
    public:
        static Matrix loadFromFile(const std::string& filepath, size_t ncols, char delimiter='\t'); 
};