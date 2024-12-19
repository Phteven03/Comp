#include "TridiagonalMatrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

TridiagonalMatrix::TridiagonalMatrix(size_t size) {

    //matrix Values
    mid.resize(size-2);
    upper.resize(size - 3);
    rightVector.resize(size - 2);

    //LU Values
    LLower.resize(size-3);
    UMid.resize(size-2);
    UUpper.resize(size-3);
}

void TridiagonalMatrix::setValues_(const std::vector<double>& x, const std::vector<double>& y) {

    //matrix Values
    size_t n = x.size();
    std::vector<double> b(n - 1);
    std::vector<double> h(n - 1);

    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
        b[i] = (y[i + 1] - y[i]) / h[i];
        if (i < n - 2) {
            upper[i] = h[i];
        }
    }

    for (size_t i = 0; i < n - 2; ++i) {
        mid[i] = 2 * (h[i] + h[i + 1]);
        rightVector[i] = 6 * (b[i + 1] - b[i]);
    }

    //LUDecomposition --> LU Values
    size_t k = mid.size();
    UMid[0] = mid[0];

    for (size_t i = 0; i < k - 1; ++i) {
        UUpper[i] = upper[i];
    }

    for (size_t i = 1; i < k; ++i) {
        LLower[i - 1] = upper[i - 1] / UMid[i - 1];
        UMid[i] = mid[i] - upper[i - 1] * LLower[i - 1];
    }
}

std::vector<double> TridiagonalMatrix::solveLU_(std::vector<double>& rightVector) const{
    size_t n = rightVector.size();
    std::vector<double> y(n);
    std::vector<double> z(n);

    y[0] = rightVector[0];
    for (size_t i = 1; i < n; ++i) {
        y[i] = rightVector[i] - LLower[i - 1] * y[i - 1];
    }

    z[n - 1] = y[n - 1] / UMid[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        z[i] = (y[i] - UUpper[i] * z[i + 1]) / UMid[i];
    }

    z.insert(z.begin(), 0.0);
    z.push_back(0.0);
    return z;
}

std::pair<size_t, std::vector<double>> TridiagonalMatrix::solveSOR_(std::vector<double>& rightVector, double omega, double tolerance) const{
    size_t n = rightVector.size();
    size_t N = mid.size();
    std::pair<size_t, std::vector<double>> itSol;
    if (n != N) {
        throw std::invalid_argument("rightVector wrong size!");
    }
    std::vector<double> guessSolution(n+1, 0.0);
    std::vector<double> previousSolution(n, 0.0);
    double diff = 0;
    double maxPrevious = 1;
    double maxGuess = 0;
    size_t iteration = 1;

    while (true) {

        for (size_t i = 1; i < n; ++i) {
            guessSolution[i] = (1 - omega) * previousSolution[i] + (omega / mid[i]) * (rightVector[i] - upper[i - 1] * guessSolution[i - 1] - upper[i] * previousSolution[i + 1]);
        }

        maxPrevious = *std::max_element(previousSolution.begin(), previousSolution.end());
        maxGuess = *std::max_element(guessSolution.begin(), guessSolution.end());
        double diff = std::abs(maxGuess - maxPrevious);

        if (diff <= tolerance) {
            previousSolution = guessSolution;
            ++iteration;
            break;
        }

        previousSolution = guessSolution;
        ++iteration;
    }
    itSol.first = iteration;
    itSol.second = guessSolution;
    guessSolution.insert(guessSolution.begin(), 0.0);
    guessSolution.push_back(0.0);
    return itSol;
}

