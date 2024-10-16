#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "mathfunc.h"

double factorial_(double x) {
    double result = 1.0;
    for (int i = 1; i <= x; ++i) {
        result *= i;
    }
    return result;
}

double logFactorial_(int n) {
    double logFact = 0.0;
    for (int i = 2; i <= n; ++i) {
        logFact += std::log(i);
    }
    return logFact;
}

double evaluatePoly_(std::vector<double> poly, double x) {
    double polyValue = 0.0;
    for (int i = 0; i < poly.size(); ++i) {
        polyValue += poly[i] * pow(x, i);
    }
    return polyValue;
}

std::vector<double> legendrePn_(int l) {
    std::vector<double> polfac(l + 1);

    for (int i = 0; i <= floor(l / 2.0); ++i) {
        int index = l - 2 * i;

        double logNumerator = logFactorial_(2 * l - 2 * i);
        double logDenominator = logFactorial_(l - i) + logFactorial_(l - 2 * i) + logFactorial_(i) + (l * std::log(2));

        polfac[index] = pow(-1, i) * exp(logNumerator - logDenominator);
    }

    return polfac;
}

std::vector<double> polyDiff_(std::vector<double> polyFunction) {
    std::vector<double> polyDiff(polyFunction.size() - 1);
    for (int i = 1; i < polyFunction.size(); ++i) {
        polyDiff[i - 1.0] = i * polyFunction[i];
    }
    return polyDiff;
}

std::vector<double> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit) {
    const double stepwidth = 1e-4;
    std::vector<double> significantValues;

    double previousValue = evaluatePoly_(poly, leftlimit);
    double previousX = leftlimit;

    for (double currentX = previousX + stepwidth; currentX < rightlimit; currentX += stepwidth) {
        double currentValue = evaluatePoly_(poly, currentX);
            if ((previousValue * currentValue) < 0.0) {
                double midpoint = (previousX + currentX) / 2.0;
                significantValues.push_back(midpoint);
            }
        previousValue = currentValue;
        previousX = currentX;
    }
    return significantValues;
}

std::vector<double> newtonRaphson_(std::vector<double> polyFunction, std::vector<double> guesses, int n) {

    std::vector<double> roots;
    std::vector<double> f = polyFunction;
    std::vector<double> fPrime = polyDiff_(f);
    double tolerance = 1e-6;
    int maxIterations = n;

    for (int j = 0; j < guesses.size(); ++j) {

        double guess = guesses[j];

        for (int i = 0; i < maxIterations; ++i) {

            double xNew = guess - evaluatePoly_(f, guess) / evaluatePoly_(fPrime, guess);

            if (abs(xNew - guess) < tolerance) {
                roots.push_back(xNew);
                break;
            }
            guess = xNew;
        }
    }
    return roots;
}

std::vector<double> polyRootFinder_(std::vector<double> poly, double n, double leftlimit, double rightlimit) {

    std::vector<double> guess = bracketing_(poly, leftlimit, rightlimit);
    size_t numberOfRoots = guess.size(); 
    std::vector<double> polyRoot(numberOfRoots);
    polyRoot = newtonRaphson_(poly, guess, n);
    return polyRoot;
}

std::vector<double> gaussLegendreWeight_(std::vector<double> legendrePoly, std::vector<double> roots) {
    size_t order = roots.size();
    std::vector<double> weights(order);
    for (size_t i = 0; i < order; ++i) {
        double PnDiff = evaluatePoly_(polyDiff_(legendrePoly), roots[i]);
        weights[i] = 2 / ((1 - roots[i] * roots[i]) * (PnDiff * PnDiff));
    }
    return weights;
}