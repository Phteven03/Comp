#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "mathfunc.h"


using namespace std;

// Function for the original integrand
// This function computes the value of the integrand 1 / sqrt(1 - y^4) for a given y
static double integrand_(double y) {
    return 1 / sqrt(1 - pow(y, 4));
}

// Implementation of the trapezoidal rule for numerical integration
// a and b are the integration bounds, n is the number of subdivisions, and integrand is the function to be integrated
static double trapezoidRule_(double a, double b, int n, function<double(double)> integrand) {
    double stepWidth = (b - a) / n;  // Step width for subdivision
    double trapIntegralValue = 0;  // Result of the trapezoidal integration
    for (int i = 1; i < n; ++i) {
        // The value is computed by summing the trapezoidal areas using the midpoint rule
        trapIntegralValue += ((integrand(a) + integrand(a + stepWidth)) / 2) * stepWidth;
        a += stepWidth;  // Increment the integration bound step by step
    }
    return trapIntegralValue;  // Return the result of the integration
}

// Implementation of Simpson's rule for numerical integration
// a and b are the integration bounds, n is the number of subdivisions, and integrand is the function to be integrated
static double simpsonRule_(double a, double b, int n, function<double(double)> integrand) {
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
    simpIntegralValue = (stepWidth / 3) * (integrand(a) + 4 * oddSumValue + 2 * evenSumValue + integrand(b - pow(b, (-10))));
    return simpIntegralValue;  // Return the result of the integration
}


static double gaussianQudrature_(double a, double b, double n, function<double(double)> integrand) {

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
    double a = 0;
    double b = 1;
    int n = 38;

    double gaussIntValue = gaussianQudrature_(a, b, n, integrand_);
    cout << "Gauss: " << gaussIntValue << endl;


    double trapResult = trapezoidRule_(a, b, n, integrand_);
    double simpResult = simpsonRule_(a, b, n, integrand_);

    cout << "Trap: " << trapResult << endl;
    cout << "Simp: " << simpResult << endl;

    //double transTrapResult = trapezoidRule_(-1, 1, n, [=](double t) { return transformedIntegrand_(a, b, t); });
    //double transSimpResult = simpsonRule_(-1, 1, n, [=](double t) { return transformedIntegrand_(a, b, t); });

    
    //cout <<  transSimpResult << endl;

    //double P3 = legendrePn_(0, 3);
    //cout << P3 << endl;
    //double dP3 = legendrePnDeriv_(0, 3);

    //cout << dP3 << endl;


}