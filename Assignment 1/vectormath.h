#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

template <typename T>
std::vector<T> vectorAddition_(std::vector<T> vector1, std::vector<T> vector2) {

	int vectorSize = vector1.size();
	for (int i = 0; i < vectorSize; ++i) {
		vector1[i] += vector2[i];
	}
	return vector1;
} 

template <typename T> 
std::vector<T> operator+(std::vector<T> vector1, std::vector<T> vector2) {
	return vectorAddition_(vector1, vector2);
}

template <typename T>
std::vector<T> vectorSubtraction_(std::vector<T> vector1, std::vector<T> vector2) {

	int vectorSize = vector1.size();
	for (int i = 0; i < vectorSize; ++i) {
		vector1[i] -= vector2[i];
	}
	return vector1;
}

template <typename T>
std::vector<T> operator-(std::vector<T> vector1, std::vector<T> vector2) {
	return vectorSubtraction_(vector1, vector2);
}

template <typename T>
T scalarProduct_(const std::vector<T>& vector1, const std::vector<T>& vector2) {
	size_t vectorSize = vector1.size();
	T sum = 0;
	for (size_t i = 0; i < vectorSize; ++i) {
		sum += vector1[i] * vector2[i];
	}
	return sum;
}

template <typename T, typename S> 
std::vector<T> scalarVectorProduct_(S scalar, std::vector<T> vector) {
	size_t vectorSize = vector.size();
	for (size_t i = 0; i < vectorSize; ++i) {
		vector[i] *= scalar; 
	}
	return vector; 
}

template <typename T, typename S>
std::vector<T> operator*(S scalar , std::vector<T> vector1) {
	return scalarVectorProduct_(scalar, vector1);
}


template <typename T, typename S> 
std::vector<T> scalarVectorProduct_(std::vector<T> vector,  S scalar) {
	size_t vectorSize = vector.size();
	for (size_t i = 0; i < vectorSize; ++i) {
		vector[i] *= scalar;
	}
	return vector;
}

template <typename T, typename S>
std::vector<T> operator*(std::vector<T> vector1, S scalar) {
	return scalarVectorProduct_(vector1, scalar);
}

template <typename T, typename S> 
std::vector<T> scalarVectorAddition_(std::vector<T> vector, S scalar) {
	size_t vectorSize = vector.size();
	for (size_t i = 0; i < vectorSize; ++i) {
		vector[i] += scalar;
	}
	return vector;
}

template <typename T, typename S>
std::vector<T> operator+(std::vector<T> vector, S scalar) {

	return scalarVectorAddition_(vector, scalar);
}

template <typename T, typename S>
std::vector<T> scalarVectorAddition_(S scalar, std::vector<T> vector) {
	size_t vectorSize = vector.size();
	for (size_t i = 0; i < vectorSize; ++i) {
		vector[i] += scalar;
	}
	return vector;
}

template <typename T, typename S>
std::vector<T> operator+(S scalar, std::vector<T> vector) {

	return scalarVectorAddition_(scalar, vector);
}

template <typename T>
std::vector<std::vector<T>> vectorVector2MatrixMultiplication_(std::vector<T> vector1, std::vector<T> vector2) {
	size_t n = vector1.size();
	std::vector<std::vector<T>> matrixProduct(n, std::vector<T>(n));
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			matrixProduct[i][j] = vector1[i] * vector2[j];
		}
	}
	return matrixProduct;
}