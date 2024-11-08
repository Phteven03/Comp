#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::vector<double>> readTxt2Matrix_(const std::string& fileName) {
    std::vector<std::vector<double>> columns;
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "File: " << fileName << " doesn't exist!" << std::endl;
        return columns;
    }

    std::string line;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream stream(line);
        double value;

        while (stream >> value) {
            row.push_back(value);
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