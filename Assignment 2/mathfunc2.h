#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>

#include "matplot/matplot.h"
#include "fftw3.h"
#include "matrixMath.h"

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
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
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

template<typename T>
void print(const T scalar) {
    std::cout << scalar << std::endl;
}

template <typename T>
void printVector(const std::vector<T>& vector) {
    std::cout << '[';
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << vector[i] << ' ';
    std::cout << ']' << std::endl;
}

template <typename T>
void printMatrix(const std::vector<std::vector<T>>& matrix) {
    for (const auto row : matrix) {
        for (T value : row) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
}

std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values, StepTimer* stepTimer = nullptr);

std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double frequency, StepTimer* stepTimer = nullptr);

std::vector<std::complex<double>> FFT_(const std::vector<double>& values, int signExp, StepTimer* stepTimer = nullptr);

std::vector<double> bubbleSort_(std::vector<double> vector);

std::vector<double> maxFinder_(std::vector<double> vector);

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth);

std::vector<double> powerMethod_(std::vector<std::vector<double>>& matrix, size_t& maxIterations);

double eigenValues_(std::vector<std::vector<double>>& matrix, std::vector<double>& eigenVector);
