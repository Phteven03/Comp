#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>

#include"fileUtils.h"
#include"mathfunc2.h"

std::vector<std::vector<double>> readTxt2Matrix_(const std::string& fileName) {
    std::vector<std::vector<double>> columns;
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "File: " << fileName << " doesn't exist!" << std::endl;
        return columns;
    }

    std::string line;
    bool firstLine = true;

    std::regex delimiterRegex("[,\\s]+");

    while (std::getline(file, line)) {
        std::vector<double> row;

        std::sregex_token_iterator begin(line.begin(), line.end(), delimiterRegex, -1);
        std::sregex_token_iterator end;

        for (auto it = begin; it != end; ++it) {
            try {
                row.push_back(std::stod(it->str()));
            }
            catch (const std::invalid_argument&) {
                //std::cerr << "Warning: Non-numeric value skipped in line: " << line << std::endl;
                continue;
            }
        }

        if (firstLine) {
            columns.resize(row.size());
            firstLine = false;
        }

        for (size_t i = 0; i < row.size(); ++i) {
            columns[i].push_back(row[i]);
        }
    }

    file.close();
    return columns;
}