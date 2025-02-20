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

struct Timer
{
    std::chrono::time_point<std::chrono::steady_clock> global_start, global_end;
    std::chrono::duration<float> global_duration;

    std::chrono::time_point<std::chrono::steady_clock> local_start, local_end;
    std::chrono::duration<float> local_duration;
    std::vector<float> laps;

    void start();
    void end();
    Timer();
    ~Timer();
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

/*
Computes the Discrete Fourier Transform (DFT) of a given real-valued input signal.
@param values A vector of real-valued input samples.
@param stepTimer (Optional) A pointer to a Timer object to measure the execution time of each step.
@return A vector of complex numbers representing the frequency components of the input signal.
*/
std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values, Timer* stepTimer = nullptr);

/*
Calculates the power spectrum of a given real-valued input signal
@param values A vector of real-valued input samples.
@param frequency The sampling frequency of the input signal.
@param stepTimer (Optional) A pointer to a Timer object to measure the execution time of each step.
@return A vector of pairs, where each pair represents a frequency and its corresponding power.
 */
std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double frequency, Timer* stepTimer = nullptr);

/*
Computes the Fast Fourier Transform (FFT) of a given real-valued input signal.
@param values A vector of real-valued input samples.
@param signExp Determines the direction of the transform (+1 for forward, -1 for inverse).
@param stepTimer (Optional) A pointer to a Timer object to measure the execution time of each step.
@return A vector of complex numbers representing the frequency components of the input signal.
*/
std::vector<std::complex<double>> FFT_(const std::vector<double>& values, int signExp, Timer* stepTimer = nullptr);

/*
Computes the inverse Fast Fourier Transform (iFFT) from the frequency domain to the time domain 
@param values A vector of complex numbers representing frequency components.
@return A vector of real-valued samples representing the reconstructed time-domain signal.
*/
std::vector<double> invFFT_(std::vector<std::complex<double>> values);

/*
Sorts a given vector of real numbers using the bubble sort algorithm. 
@param vector A vector of real numbers to be sorted.
@return A sorted vector in ascending order.
*/
std::vector<double> bubbleSort_(std::vector<double> vector);

/*
Finds the maximum values in a given vector of real numbers.
@param vector A vector of real numbers.
@return A vector containing the maximum values found in the input vector.
*/
std::vector<double> maxFinder_(std::vector<double> vector);

/*
Computes cubic spline interpolation for given input points and step width. 
@param x A vector of x-coordinates.
@param y A vector of y-coordinates.
@param z A vector of z-values for three-dimensional interpolation.
@param stepWidth The step size for interpolation.
@return A splineValues object containing the computed spline coefficients.
*/
splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth);

/*
Performs the power method to estimate the dominant eigenvector of a matrix. 
@param matrix A reference to a 2D vector representing the input matrix.
@param maxIterations A reference to the maximum number of iterations to perform.
@return A vector representing the estimated dominant eigenvector.
*/
std::vector<double> powerMethod_(std::vector<std::vector<double>>& matrix, size_t& maxIterations);

/*
Computes the dominant eigenvalue of a matrix given its corresponding eigenvector. 
@param matrix A reference to a 2D vector representing the input matrix.
@param eigenVector A reference to the corresponding eigenvector.
@return The computed dominant eigenvalue of the matrix.
*/
double eigenValues_(std::vector<std::vector<double>>& matrix, std::vector<double>& eigenVector);
