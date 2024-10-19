#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>
#include "vectormath.h"
#include "matplot/matplot.h"

//exercise 1
#include "mathfunc1a.h"
#include "plot1a.h"
#include "mathfunc2a.h"

double PI = 3.14159265359;

//exercise 1a

static double integrand1a_(double y) {
    return 1 / std::sqrt(1 - std::pow(y, 4));
}

static double trapezoidRule_(double a, double b, int n, std::function<double(double)> integrand) {
    double stepWidth = (b - a) / n;
    double trapIntegralValue = 0;
    for (int i = 1; i < n; ++i) {
        trapIntegralValue += ((integrand(a) + integrand(a + stepWidth)) / 2) * stepWidth;
        a += stepWidth;
    }
    return trapIntegralValue;
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
    return simpIntegralValue;
}

static double gaussianQudrature_(double a, double b, double n, std::function<double(double)> integrand) {

    double gaussIntValue = 0;

    std::vector<double> weights(n);
    std::vector<double> roots(n);


    roots = legendreRootFinder_(-1, 1, n);
    weights = gaussLegendreWeight_(roots, n);

    for (int i = 0; i < n; ++i) {
        double transform = ((b - a) / 2) * roots[i] + ((a + b) / 2);
        gaussIntValue += weights[i] * integrand(transform);
    }
    gaussIntValue = (b - a) / 2 * gaussIntValue;
    return gaussIntValue;
}

//exercise 1b

static double integrand1b1_(double x, double a) {
    return (std::sqrt(cosh(a) - std::cosh(x))) / (cosh(a) - std::cosh(x));
}
static double integrand1b2_(double x, double a) {
    return 1 / (std::exp(a) - std::exp(abs(x)));
}
static double integrand1b3_(double x, double a) {
    return 1 / std::sqrt(std::cos(x) - std::cos(a));
}

int main() {

    //exercise 1a
    const double a = 0;
    const double b = 1;
    const double n = 100;
    const double REALINTVALUE = 1.311028777146120;

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
        gaussIntValue = gaussianQudrature_(a, b, i, integrand1a_);
        gaussIntVector.push_back(gaussIntValue);

        trapIntValue = trapezoidRule_(a, b, i, integrand1a_);
        trapIntVector.push_back(trapIntValue);

        simpIntValue = simpsonRule_(a, b, i, integrand1a_);
        simpIntVector.push_back(simpIntValue);
    }

    //plotresult1a_(trapIntVector, simpIntVector, gaussIntVector, n);


    //exercise 1b
    std::vector<std::vector<double>> allIntegrals;
    for (double a = 0; a <= PI; a += (PI / 1e2)) {

        std::vector<std::function<double(double)>> integrandFixedA = {
            [a](double x) { return integrand1b1_(x, a); },
            [a](double x) { return integrand1b2_(x, a); },
            [a](double x) { return integrand1b3_(x, a); }
        };

        std::vector<double> intVectorForA;

        for (const auto& integrand : integrandFixedA) {
            double intVal = gaussianQudrature_(0, a, 50, integrand);
            intVectorForA.push_back(intVal);
        }

        allIntegrals.push_back(intVectorForA);
    }

    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b1_: " << allIntegrals[i][0] << std::endl;
    }
    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b2_: " << allIntegrals[i][1] << std::endl;
    }
    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b3_: " << allIntegrals[i][2] << std::endl;
    }


   

    /*const int n = 100;
    plotresult1a_(n);*/



}