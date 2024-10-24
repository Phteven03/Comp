#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "mathfunc1a.h"
#include "vectormath.h"
#include <chrono>

// Timer constructor: starts the timer
Timer::Timer() {
    start = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(0.0f);
    end = start;
}

// Timer destructor: ends the timer and prints the duration in milliseconds
Timer::~Timer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    float ms = duration.count() * 1000.0f;
    std::cout << "Timer took " << ms << "ms " << std::endl;
}

// Computes Legendre polynomial P_n(x) using recursion
double legendrePn_(double x, int n) {
    double Pn_1 = 1;
    double Pn = x;

    if (n == 0) {
        return Pn_1; // Return P_0(x)
    }
    else if (n == 1) {
        return Pn; // Return P_1(x)
    }

    // Recursive calculation for higher degrees
    for (int i = 2; i <= n; ++i) {
        double Pn_next = ((2 * i - 1) * x * Pn - (i - 1) * Pn_1) / i;
        Pn_1 = Pn;
        Pn = Pn_next;
    }

    return Pn;
}

// Computes the derivative of Legendre polynomial P_n'(x)
double legendrePnDiff_(double x, int n) {
    if (n == 0) {
        return 0; // Derivative of P_0(x) is 0
    }
    else if (n == 1) {
        return 1; // Derivative of P_1(x) is 1
    }
    else {
        // Derivative formula for higher degrees
        return (n / (1 - x * x)) * (legendrePn_(x, n - 1) - x * legendrePn_(x, n));
    }
}

// Uses bracketing to find intervals where roots of P_n(x) may exist
std::vector<double> legendreBracketing_(double leftlimit, double rightlimit, double n) {
    std::vector<double> guesses;
    const double stepwidth = 1e-5;

    double previousValue = legendrePn_(leftlimit, n);
    double previousX = leftlimit;

    // Loop through the interval and check for sign changes
    for (double currentX = previousX + stepwidth; currentX < rightlimit; currentX += stepwidth) {
        double currentValue = legendrePn_(currentX, n);
        if ((previousValue * currentValue) < 0.0) { // Sign change indicates a root
            double midpoint = (previousX + currentX) / 2.0;
            guesses.push_back(midpoint); // Store midpoint as a guess for root
        }
        previousValue = currentValue;
        previousX = currentX;
    }
    return guesses;
}

// Uses Newton-Raphson method to find roots of P_n(x)
std::vector<double> legendreNewtonRaphson_(std::vector<double> guesses, double n) {
    std::vector<double> roots;

    double tolerance = 1e-6;
    int maxIterations = 100;

    // Loop over each initial guess to find the root
    for (size_t j = 0; j < guesses.size(); ++j) {
        double guess = guesses[j];
        for (int i = 0; i < maxIterations; ++i) {
            double f = legendrePn_(guess, n);
            double fPrime = legendrePnDiff_(guess, n);

            double xNew = guess - f / fPrime; // Newton-Raphson update

            if (abs(xNew - guess) < tolerance) { // Convergence check
                roots.push_back(xNew); // Store the converged root
                break;
            }
            guess = xNew;
        }
    }
    return roots;
}

// Finds the roots of the Legendre polynomial in the given interval
std::vector<double> legendreRootFinder_(double leftlimit, double rightlimit, double n) {
    std::vector<double> guesses = legendreBracketing_(leftlimit, rightlimit, n);
    std::vector<double> roots = legendreNewtonRaphson_(guesses, n);
    return roots;
}

// Computes Gauss-Legendre quadrature weights based on the roots of P_n(x)
std::vector<double> gaussLegendreWeight_(std::vector<double> roots, double n) {
    size_t order = roots.size();
    std::vector<double> weights(order);

    // Calculate weights for each root
    for (size_t i = 0; i < order; ++i) {
        double PnDiff = legendrePnDiff_(roots[i], n);
        weights[i] = 2 / ((1 - roots[i] * roots[i]) * (PnDiff * PnDiff)); // Formula for Gauss-Legendre weight
    }
    return weights;
}
