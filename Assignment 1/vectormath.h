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
std::vector<T> scalarProduct_(std::vector<T> vector1, std::vector<T> vector2) {
	
	size_t vectorSize = vector1.size();
	for (size_t i = 0; i < vectorSize; ++i) {
		vector1[i] *= vector2[i];
	}
	T sum = accumulate(vector1.begin(), vector1.end(), 0);
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