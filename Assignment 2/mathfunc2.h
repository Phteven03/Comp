#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>

#include "matplot/matplot.h"
#include "fftw3.h"

struct Timer {
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::duration<float> duration;

    Timer();

    ~Timer();
};

class StepTimer {
public:
    void startTimer();
    void stopStoreTimer();
    std::vector<float> getTimes();
private:
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<float> duration;
    std::vector<float> times;
};


struct TridiagonalMatrix {
    std::vector<double> mid;
    std::vector<double> upper;
    std::vector<double> rightvector;
};

struct LUDecomposition {
    std::vector<double> LLower;
    std::vector<double> UMid;
    std::vector<double> UUpper;
};

struct splineValues {
    std::vector<double> xValues;
    std::vector<double> splineValues;
};

struct SOR {
    std::vector<double> iterations;
    std::vector<double> solution;

};



template <typename T>
void printVector(const std::vector<T>& vector) {
    std::cout << '[';
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << vector[i] << ' ';
    std::cout << ']' << std::endl;
}

std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values, StepTimer* stepTimer = nullptr);

std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double frequency, StepTimer* stepTimer = nullptr);

std::vector<std::complex<double>> FFT_(const std::vector<double>& values, StepTimer* stepTimer = nullptr);
 
TridiagonalMatrix createTridiagonalMatrix_(std::vector<double>& x, std::vector<double>& y);

LUDecomposition LUD_(TridiagonalMatrix& matrix);

std::vector<double> solveLU_(LUDecomposition& LU, std::vector<double>& u);

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth);

SOR calculateWithSOR_(TridiagonalMatrix& matrix, double& omega, double tolerance);