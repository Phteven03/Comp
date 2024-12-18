#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

class TridiagonalMatrix {
public:
    std::vector<double> mid;
    std::vector<double> upper;
    std::vector<double> rightvector;
    std::vector<double> LLower;
    std::vector<double> UMid;
    std::vector<double> UUpper;
public:
    TridiagonalMatrix(size_t size);
    void setValues(const std::vector<double>& x, const std::vector<double>& y);
    TridiagonalMatrix LUD_() const;
    std::vector<double> solveLU_() const;

};