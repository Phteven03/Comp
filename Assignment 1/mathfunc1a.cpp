#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "mathfunc1a.h"
#include "vectormath.h"
#include <chrono>

Timer::Timer() {
    start = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(0.0f);
    end = start;
}

Timer::~Timer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    float ms = duration.count() * 1000.0f;
    std::cout << "Timer took " << ms << "ms " << std::endl;
}

double legendrePn_(double x, int n) {
    double Pn_1 = 1;
    double Pn = x;

    if (n == 0) {
        return Pn_1;
    }
    else if (n == 1) {
        return Pn;
    }

    for (int i = 2; i <= n; ++i) {
        double Pn_next = ((2 * i - 1) * x * Pn - (i - 1) * Pn_1) / i;
        Pn_1 = Pn;
        Pn = Pn_next;
    }

    return Pn;
}

double legendrePnDiff_(double x, int n) {
    if (n == 0) {
        return 0;
    }
    else if (n == 1) {
        return 1;
    }
    else {
        return (n / (1 - x * x)) * (legendrePn_(x, n - 1) - x * legendrePn_(x, n));
    }
}

std::vector<double> legendreBracketing_(double leftlimit, double rightlimit , double n) {
    const double stepwidth = 1e-5;
    double significantValue = 0;
    std::vector<double> guesses;

    double previousValue = legendrePn_(leftlimit , n);
    double previousX = leftlimit;

    for (double currentX = previousX + stepwidth; currentX < rightlimit; currentX += stepwidth) {
        double currentValue = legendrePn_(currentX,n);
        if ((previousValue * currentValue) < 0.0) {
             double midpoint = (previousX + currentX) / 2.0;
             guesses.push_back(midpoint);
        }
        previousValue = currentValue;
        previousX = currentX;
    }
    return guesses;

}

std::vector<double> legendreNewtonRaphson_(std::vector<double> guesses, double n) {

    std::vector<double> roots;

    size_t sizeGuesses = guesses.size();
    double tolerance = 1e-6;
    int maxIterations = 100;

    for (size_t j = 0; j < sizeGuesses; ++j) {
        double guess = guesses[j];
        for (int i = 0; i < maxIterations; ++i) {

            double f = legendrePn_(guess, n);
            double fPrime = legendrePnDiff_(guess, n);

            double xNew = guess - f / fPrime;

            if (abs(xNew - guess) < tolerance) {
                roots.push_back(xNew);
                break;
            }
            else if (i == maxIterations) {
                roots.push_back(xNew);
                break;
            }
            guess = xNew;
        }
    }
    return roots;
}

std::vector<double> legendreRootFinder_(double leftlimit, double rightlimit, double n) {

    std::vector<double> guesses = legendreBracketing_(leftlimit, rightlimit , n);
    std::vector<double> roots = legendreNewtonRaphson_(guesses , n);
    return roots;
}

std::vector<double> gaussLegendreWeight_(std::vector<double> roots, double n) {
    size_t order = roots.size();
    std::vector<double> weights(order);
    for (size_t i = 0; i < order; ++i) {
        double PnDiff = legendrePnDiff_(roots[i],n);
        weights[i] = 2 / ((1 - roots[i] * roots[i]) * (PnDiff * PnDiff));
    }
    return weights;
} 