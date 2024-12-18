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

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth);
