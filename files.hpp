#pragma once

#include<filesystem>

class MFile {
    private:
        std::filesystem::path mPath{};
        std::map<double, std::complex<double>> mData{};
        const char mDelimiter = '\t';
    public:
        explicit MFile(const std::filesystem::path& path);
        MFile(const std::filesystem::path& path, const char delimiter);

        int load();
        std::map<double, std::complex<double>> getData() const;
};