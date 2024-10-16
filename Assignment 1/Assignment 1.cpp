#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "mathfunc.h"
#include "matplot/matplot.h"
#include "plots.h"
#include "vectormath.h"

// Function for the original integrand
// This function computes the value of the integrand 1 / sqrt(1 - y^4) for a given y
static double integrand_(double y) {
    return 1 / std::sqrt(1 - std::pow(y, 4));
}

// Implementation of the trapezoidal rule for numerical integration
// a and b are the integration bounds, n is the number of subdivisions, and integrand is the function to be integrated
static double trapezoidRule_(double a, double b, int n, std::function<double(double)> integrand) {
    double stepWidth = (b - a) / n;  // Step width for subdivision
    double trapIntegralValue = 0;  // Result of the trapezoidal integration
    for (int i = 1; i < n; ++i) {
        // The value is computed by summing the trapezoidal areas using the midpoint rule
        trapIntegralValue += ((integrand(a) + integrand(a + stepWidth)) / 2) * stepWidth;
        a += stepWidth;  // Increment the integration bound step by step
    }
    return trapIntegralValue;  // Return the result of the integration
}

static double simpsonRule_(double a, double b, int n, std::function<double(double)> integrand) {
    double simpIntegralValue = 0;  // Result of Simpson's integration
    double stepWidth = (b - a) / n;  // Step width for subdivision
    double oddSumValue = 0;  // Sum of function values at odd indices
    double evenSumValue = 0;  // Sum of function values at even indices

    // Loop over odd indices to sum function values
    for (int i = 1; i < n; i += 2) {
        oddSumValue += integrand(a + i * stepWidth);
    }

    // Loop over even indices to sum function values
    for (int i = 2; i < n; i += 2) {
        evenSumValue += integrand(a + i * stepWidth);
    }

    // Simpson's formula, which combines the function values at the boundaries, odd, and even points
    simpIntegralValue = (stepWidth / 3) * (integrand(a) + 4 * oddSumValue + 2 * evenSumValue + integrand(b - std::pow(b, (-10))));
    return simpIntegralValue;  // Return the result of the integration
}

static double gaussianQudrature_(double a, double b, double n, std::function<double(double)> integrand) {

    double gaussIntValue = 0;
    std::vector<double> weights(n);
    std::vector<double> roots(n);
    std::vector<double> legPoly(n + 1);

    legPoly = legendrePn_(n);
    roots = polyRootFinder_(legPoly, n, -1, 1);
    weights = gaussLegendreWeight_(legPoly,roots);


    for (int i = 0; i < n; ++i) {
        double transform = ((b - a) / 2) * roots[i] + ((a + b) / 2);
        gaussIntValue += weights[i] * integrand(transform);
    }
    gaussIntValue = (b - a) / 2 * gaussIntValue;
    return gaussIntValue;
}


int main() {
    const double a = 0;
    const double b = 1;
    const int n = 30;
    const double REALINTVALUE = 1.311028777146120;


    std::vector<int> nVector;
    for (int i = 1; i <= n; ++i) {
        nVector.push_back(i);
    }

    std::vector<double> gaussIntVector;
    std::vector<double> trapIntVector;
    std::vector<double> simpIntVector;

    std::vector<double> gaussIntError;
    std::vector<double> trapIntError;
    std::vector<double> simpIntError;

    double trapIntValue = 0.0;
    double gaussIntValue = 0.0;
    double simpIntValue = 0.0;

    double gaussIntErrorV = 0.0;
    double trapIntErrorV = 0.0;
    double simpIntErrorV = 0.0;

    for (int i = 0; i < n; ++i) {
        gaussIntValue = gaussianQudrature_(a, b, i, integrand_);
        gaussIntVector.push_back(gaussIntValue);

        trapIntValue = trapezoidRule_(a, b, i, integrand_);
        trapIntVector.push_back(trapIntValue);

        simpIntValue = simpsonRule_(a, b, i, integrand_);
        simpIntVector.push_back(simpIntValue);
    }

    matplot::plot(nVector, gaussIntVector);
    matplot::hold(matplot::on);
    matplot::plot(nVector, trapIntVector);
    matplot::plot(nVector, simpIntVector);
    matplot::xlabel("N");
    matplot::ylabel("Error");
    matplot::show();
}