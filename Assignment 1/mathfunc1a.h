#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>

/*
  Timer struct for measuring execution time.
 */
struct Timer {
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::duration<float> duration;

    Timer();

    ~Timer();
};

/*
  Function to print the contents of a vector.
 
  @param vector: A reference to a vector of type T to be printed.
 */
template <typename T>
void printVector(std::vector<T>& vector) {
    std::cout << '[';
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << vector[i] << ' ';
    std::cout << ']' << std::endl;
}

/*
  Evaluates the Legendre polynomial of degree n at a given point x.
 
  @param x: The point at which the polynomial is evaluated.
  @param n: The degree of the Legendre polynomial.
  @return The value of the Legendre polynomial P_n(x).
 */
double legendrePn_(double x, int n);

/*
  Evaluates the derivative of the Legendre polynomial of degree n at a given point x.
 
  @param x: The point at which the derivative is evaluated.
  @param n: The degree of the Legendre polynomial.
  @return The derivative of the Legendre polynomial P_n'(x).
 */
double legendrePnDiff_(double x, int n);

/*
  Finds the roots of the Legendre polynomial using the Newton-Raphson method.
 
  @param guesses: A vector of initial guesses for the roots.
  @param n: The degree of the Legendre polynomial.
  @return A vector of converged root values.
 */
std::vector<double> legendreNewtonRaphson_(std::vector<double> guesses, double n);

/*
  Uses bracketing to find intervals where the roots of the Legendre polynomial may exist.
 
  @param leftlimit: The left bound of the interval.
  @param rightlimit: The right bound of the interval.
  @param n: The degree of the Legendre polynomial.
  @return A vector of intervals containing potential roots.
 */
std::vector<double> legendreBracketing_(double leftlimit, double rightlimit, double n);

/*
  Computes the weights for Gauss-Legendre quadrature using the roots of the Legendre polynomial.
 
  @param roots: A vector containing the roots of the Legendre polynomial.
  @param n: The degree of the Legendre polynomial.
  @return A vector of quadrature weights corresponding to the roots.
 */
std::vector<double> gaussLegendreWeight_(std::vector<double> roots, double n);

/*
  Finds the roots of the Legendre polynomial in a given interval using a combination of bracketing and root-finding techniques.
 
  @param leftlimit: The left bound of the interval.
  @param rightlimit: The right bound of the interval.
  @param n: The degree of the Legendre polynomial.
  @return A vector of the roots of the Legendre polynomial within the specified interval.
 */
std::vector<double> legendreRootFinder_(double leftlimit, double rightlimit, double n);
