#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <vectormath.h>

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

struct charge {
    double theta;
    double phi;
    double gradPhi;
    double gradTheta;
};

double orthodrome_(const charge& q1, const charge& q2);

double totalEnergy_(std::vector<charge> charges);

void numericalGradient_(std::vector<charge>& charges, double stepSize);

void gradientDecent_(std::vector<charge>& charges, double stepSize, double maxIterations);


std::vector<long double> gravitationalForce_(std::vector<long double>& r, std::vector<long double>& R, long double massLight, long double massHeavy);

std::vector<long double> inertialForce_(std::vector<long double>& r, std::vector<long double>& velocity, std::vector<long double>& omega, long double massLight);

std::vector<long double> totalForce_(std::vector<long double>& r, std::vector<long double>& R1, std::vector<long double>& R2, std::vector<long double>& velocity, std::vector<long double>& omega, long double m, long double M1, long double M2);

std::vector<long double> totalForceDimLess_(std::vector<long double>& r, std::vector<long double>& R1, std::vector<long double>& R2, std::vector<long double>& velocity, std::vector<long double>& omega, long double m, long double M1, long double M2);

std::vector<std::vector<long double>> forwardEuler_(std::vector<long double> r, std::vector<long double>& R1, std::vector<long double>& R2, std::vector<long double> velocity, std::vector<long double>& omega, long double m, long double M1, long double M2, long double dt, int n);


