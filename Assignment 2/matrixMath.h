#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

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
    std::vector<T> matrixVectorProduct(n, 0);
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

template <typename T>
std::vector<T> cumsumVector(std::vector<T> vector)
{
    size_t n = vector.size();
    T sum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += vector[i];
        vector[i] = sum;
    }

    return vector;
}