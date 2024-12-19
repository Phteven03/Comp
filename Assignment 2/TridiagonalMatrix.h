#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

class TridiagonalMatrix {
private:
    std::vector<double> mid;
    std::vector<double> upper;
    std::vector<double> LLower;
    std::vector<double> UMid;
    std::vector<double> UUpper;
public:
    std::vector<double> rightVector;
    TridiagonalMatrix(size_t size);
    void setValues_(const std::vector<double>& x, const std::vector<double>& y);
    std::vector<double> solveLU_(std::vector<double>& rightVector) const;
    std::pair<size_t, std::vector<double>> solveSOR_(std::vector<double>& rightVector, double omega, double tolerance) const;
};