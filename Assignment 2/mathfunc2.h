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

template <typename T>
void printMatrix(const std::vector<std::vector<T>>& matrix) {
    for(const auto row : matrix) {
        for (T value : row) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
}

std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values, StepTimer* stepTimer = nullptr);

std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double frequency, StepTimer* stepTimer = nullptr);

std::vector<std::complex<double>> FFT_(const std::vector<double>& values, StepTimer* stepTimer = nullptr);

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth);

template <typename T>
std::vector<std::vector<T>> matrixMultiplication_(std::vector<std::vector<T>> matrixA, std::vector<std::vector<T>> matrixB) {
    size_t n = matrixA[0].size();
    std::vector<std::vector<T>> matrixProduct(n, std::vector<T>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < n; ++k) {
            for (size_t j = 0; j < n; ++j) {
                matrixProduct[i][k] += matrixA[i][j] * matrixB[j][k];
            }
        }
    }
    return matrixProduct;
}

template <typename T> 
std::vector<T> matrixVectorMuliplicaton_(std::vector<std::vector<T>> matrix, std::vector<T> vector) {
    size_t n = vector.size();
    std::vector<T> matrixVectorProduct(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrixVectorProduct[i] += matrix[i][j] * vector[j];
        }
    }
    return matrixVectorProduct;
}

template <typename T>
std::vector<std::vector<T>> scalarMatrixMultiplication_(T scalar, std::vector<std::vector<T>> matrix) {
    size_t n = matrix[0].size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix[i][j] *= scalar;
        }
    }
    return matrix;
}

template <typename T>
std::vector<std::vector<T>> matrixMatrixAddition_(std::vector<std::vector<T>> matrix1, std::vector<std::vector<T>> matrix2) {
    size_t n = matrix1[1].size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix1[i][j] += matrix2[i][j];
        }
    }
    return matrix1;
}

template <typename T>
std::vector<std::vector<T>> matrixMatrixSubtraction_(std::vector<std::vector<T>> matrix1, std::vector<std::vector<T>> matrix2) {
    size_t n = matrix1[1].size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix1[i][j] -= matrix2[i][j];
        }
    }
    return matrix1;
}

std::vector<double> powerMethod_(std::vector<std::vector<double>>& matrix, size_t& maxIterations);

double eigenValues_(std::vector<std::vector<double>>& matrix, std::vector<double>& eigenVector);
