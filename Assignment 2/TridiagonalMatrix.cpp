#include "TridiagonalMatrix.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

TridiagonalMatrix::TridiagonalMatrix(size_t size) {
    mid.resize(size-2);
    upper.resize(size - 3);
    rightvector.resize(size-2);
    LLower.resize(size-3);
    UMid.resize(size-2);
    UUpper.resize(size-3);
}

void TridiagonalMatrix::setValues(const std::vector<double>& x, const std::vector<double>& y) {
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
        rightvector[i] = 6 * (b[i + 1] - b[i]);
    }
}

TridiagonalMatrix TridiagonalMatrix::LUD_() const {
    TridiagonalMatrix matrix = *this;
    size_t n = mid.size();

    matrix.UMid[0] = matrix.mid[0];

    for (size_t i = 0; i < n - 1; ++i) {
        matrix.UUpper[i] = matrix.upper[i];
    }

    for (size_t i = 1; i < n; ++i) {
        matrix.LLower[i - 1] = matrix.upper[i - 1] / matrix.UMid[i - 1];
        matrix.UMid[i] = matrix.mid[i] - matrix.upper[i - 1] * matrix.LLower[i - 1];
    }
    return matrix;
}

std::vector<double> TridiagonalMatrix::solveLU_() const{
    TridiagonalMatrix matrix = *this;
    size_t n = matrix.rightvector.size();
    std::vector<double> y(n);
    std::vector<double> z(n);

    y[0] = matrix.rightvector[0];
    for (size_t i = 1; i < n; ++i) {
        y[i] = matrix.rightvector[i] - matrix.LLower[i - 1] * y[i - 1];
    }

    z[n - 1] = y[n - 1] / matrix.UMid[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        z[i] = (y[i] - matrix.UUpper[i] * z[i + 1]) / matrix.UMid[i];
    }

    z.insert(z.begin(), 0.0);
    z.push_back(0.0);
    return z;
}