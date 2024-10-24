#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>
#include "vectormath.h"
#include "matplot/matplot.h"

//exercise 1
#include "mathfunc1a.h"
#include "plot1.h"

//exercise 2
#include "mathfunc2a.h"

//exercise 3
#include <random>

const double PI = 3.14159265359;
// Small error term to avoid singularities
const double SINGULARITY_ERROR_TERM = 1e-5; 


// ---------    exercise 1a   -----------



// Integrand function for 1a
static double integrand1a_(double y) {
    return 1 / std::sqrt(1 - std::pow(y, 4));
}

// Trapezoid rule implementation for numerical integration
static double trapezoidRule_(double a, double b, int n, std::function<double(double)> integrand) {
    if (isinf(integrand(a))) a += SINGULARITY_ERROR_TERM; // Handle singularities at boundaries
    if (isinf(integrand(b))) b -= SINGULARITY_ERROR_TERM;
    double stepWidth = (b - a) / n;
    double trapIntegralValue = 0.5 * (integrand(a) + integrand(b));
    for (int i = 1; i < n; ++i) {
        trapIntegralValue += integrand(a + i * stepWidth);
    }
    return trapIntegralValue * stepWidth;
}

// Simpson's rule implementation for numerical integration
static double simpsonRule_(double a, double b, int n, std::function<double(double)> integrand) {
    double simpIntegralValue = 0;
    double stepWidth = (b - a) / n;
    double oddSumValue = 0;
    double evenSumValue = 0;

    for (int i = 1; i < n; i += 2) {
        oddSumValue += integrand(a + i * stepWidth); // Odd indexed terms
    }
    for (int i = 2; i <= n; i += 2) {
        evenSumValue += integrand(a + i * stepWidth); // Even indexed terms
    }

    // Apply Simpson's rule formula
    simpIntegralValue = (stepWidth / 3) * (integrand(a) + 4 * oddSumValue + 2 * evenSumValue + integrand(b - SINGULARITY_ERROR_TERM));
    return simpIntegralValue;
}

// Gaussian quadrature implementation using Legendre polynomials
static double gaussianQudrature_(double a, double b, double n, std::function<double(double)> integrand) {
    double gaussIntValue = 0;
    std::vector<double> weights(n);
    std::vector<double> roots(n);

    // Find roots and weights of Legendre polynomials
    roots = legendreRootFinder_(-1, 1, n);
    weights = gaussLegendreWeight_(roots, n);

    // Apply Gaussian quadrature formula
    for (int i = 0; i < n; ++i) {
        double transform = ((b - a) / 2) * roots[i] + ((a + b) / 2); // Transform roots to match interval [a, b]
        gaussIntValue += weights[i] * integrand(transform);
    }
    gaussIntValue = (b - a) / 2 * gaussIntValue;
    return gaussIntValue;
}




// ---------    exercise 1b   -----------

// Three integrand functions for 1b, parameterized by 'a'
static double integrand1b1_(double x, double a) {
    return (std::sqrt(cosh(a) - std::cosh(x))) / (cosh(a) - std::cosh(x));
}

static double integrand1b2_(double x, double a) {
    return 1 / std::sqrt((std::exp(a) - std::exp(x)));
}

static double integrand1b3_(double x, double a) {
    return 1 / std::sqrt(std::cos(x) - std::cos(a));
}



// ---------    exercise 1c   -----------

// Potential function for 1c
static double potential_(double k, double x, double d) {
    return std::tanh(k * x);
}

// Two integrand functions for 1c, using different scalings of 'k'
static double integrand1c1_(double x, double k) {
    return 1 / std::sqrt(std::abs(std::tanh(k)) - std::abs(std::tanh(k * x)));
}

static double integrand1c2_(double x, double k) {
    return 1 / std::sqrt(std::abs(std::tanh(0.5 * k)) - std::abs(std::tanh(k * x)));
}



// ---------    exercise 2a   -----------

double R = 14959787000; //Distance earth sun
double mu = 3e-6; // \frac{M_2}{M_1+M_2}
// Coefficients for polynomial in exercise 2a
std::vector<double> function2a = { -mu, 2 * mu, -mu, 3 - 2 * mu, mu - 3, 1 }; 





int main() {

    // ---------    exercise 1a   -----------
    /*
        
    const double a = 0;
    const double b = 1;
    const double n = 100;

    std::vector<double> gaussIntVector, trapIntVector, simpIntVector;
    std::vector<double> gaussIntError, trapIntError, simpIntError;

    double trapIntValue = 0.0, gaussIntValue = 0.0, simpIntValue = 0.0;
    double gaussIntErrorV = 0.0, trapIntErrorV = 0.0, simpIntErrorV = 0.0;

    // Calculate integrals using different methods
    for (int i = 1; i < n; ++i) {
        gaussIntValue = gaussianQudrature_(a, b, i, integrand1a_);
        gaussIntVector.push_back(gaussIntValue);

        trapIntValue = trapezoidRule_(a, b, i, integrand1a_);
        trapIntVector.push_back(trapIntValue);
    }
    for (int i = 1; i < 2 * n; i += 2) {
        simpIntValue = simpsonRule_(a, b, i, integrand1a_);
        simpIntVector.push_back(simpIntValue);
    }

    plotresult1a_(trapIntVector, simpIntVector, gaussIntVector, n); // Plot results for 1a
    */


    // ---------    exercise 1b   -----------

    /*
    std::vector<std::vector<double>> allIntegrals;
    int n1b = 50;

    // Calculate integrals for 1b
    for (double a = 0; a <= PI; a += (PI / 1e2)) {
        std::vector<std::function<double(double)>> integrandFixedA = {
            [a](double x) { return integrand1b1_(x, a); },
            [a](double x) { return integrand1b2_(x, a); },
            [a](double x) { return integrand1b3_(x, a); }
        };

        std::vector<double> intVectorForA;

        for (const auto& integrand : integrandFixedA) {
            double intVal = gaussianQudrature_(0, a, n1b, integrand);
            intVectorForA.push_back(intVal);
        }

        allIntegrals.push_back(intVectorForA);
    }

    plotresult1b_(allIntegrals); // Plot results for 1b

    // Print integral values for 1b
    //for (size_t i = 0; i < allIntegrals.size(); ++i) {
    //    std::cout << "Result for a[" << i << "] with integrand1b1_: " << allIntegrals[i][0] << std::endl;
    //    std::cout << "Result for a[" << i << "] with integrand1b2_: " << allIntegrals[i][1] << std::endl;
    //    std::cout << "Result for a[" << i << "] with integrand1b3_: " << allIntegrals[i][2] << std::endl;
    //}
    */
    


    // ---------    exercise 1c   -----------
    
    /*
    std::vector<double> potential1c1Vector;
    std::vector<double> potential1c05Vector;
    std::vector<double> potentialfracVec;
    std::vector<double> kVector;

    // Calculate potentials for 1c
    for (double k = 0; k <= 3; k += 1e-2) {
        kVector.push_back(k);
        auto integrandFixedK = [k](double x) { return integrand1c1_(x, k); };
        double potential1c1 = gaussianQudrature_(0, 1, 50, integrandFixedK);
        potential1c1Vector.push_back(potential1c1);

        auto integrandFixedK05 = [k](double x) { return integrand1c2_(x, k); };
        double potential1c05 = gaussianQudrature_(0, 0.5, 50, integrandFixedK05);
        potential1c05Vector.push_back(potential1c05);

        double potentialfrac = potential1c1 / potential1c05; // Ratio of potentials
        potentialfracVec.push_back(potentialfrac);
    }

    plotresult1c_(kVector, potential1c1Vector, potential1c05Vector);
    */


    // ---------    exercise 2a   -----------
    /*
    
    //!!!!!!!!!!!!!!change stepwidth of bracketing!!!!!!!!!!!!!!!
    
    std::vector<double> convergenceNewton = newtonConvergence_(function2a, -3, 3);
    std::vector<double> convergenceBisection = bisectonConvergence_(function2a, -3, 3);
    

    std::vector<double> rootsbis = polyRootBisection_(function2a, -3, 3);
    std::vector<double> distanceL1Bi = scalarVectorProduct_(R, rootsbis);
    std::cout << "Ratio: " << " ";
    printVector(rootsbis);
    std::cout << "Bisection Distance L1 (in meters):" << " ";
    printVector(distanceL1Bi);


    std::vector<double> rootsnet = polyRootNewtonRaphson_(function2a, -3, 3);
    std::vector<double> distanceL1New = scalarVectorProduct_(R, rootsnet);
    std::cout << "Ratio: " << " ";
    printVector(rootsnet);
    std::cout << "Newton Distance L1 (in meters):" << " ";
    printVector(distanceL1New);
    


    plotresult2a_(convergenceBisection, convergenceNewton); // Plot convergence results for 2a
    */


    // ---------    exercise 2b   ----------- 
    /*
    //remove the comments from task 2a
    
    std::vector<double> errorIVec, errorIp1Vec;

    // Calculate errors for 2b
    for (size_t i = 0; i < convergenceNewton.size(); ++i) {
        double trueRoot = convergenceNewton.back();
        double errorI = std::abs(convergenceNewton[i] - trueRoot);
        errorIVec.push_back(errorI);
        double errorIp1 = std::abs(convergenceNewton[i + 1] - trueRoot);
        errorIp1Vec.push_back(errorIp1);
    }

    // Calculate and print convergence rate
    double k = (errorIVec[1] - errorIVec.back()) / (errorIp1Vec[1] - errorIp1Vec.back());
    std::cout << k << std::endl;
    plotresult2b_(errorIVec, errorIp1Vec);
    */

    // ---------    exercise 3ab   -----------
    
    /*
    //!!!!!!!!!!!!!!change stepwidth of bracketing!!!!!!!!!!!!!!!
    std::vector<std::vector<double>> rootsVec;
    std::vector<double> rootsNumberVec;

    // Generate random polynomials and find their roots
    for (int i = 0; i < 1000; ++i) {
        std::vector<double> randpoly = generateRandomNumbers(7);
        std::vector<double> roots = polyRootNewtonRaphson_(randpoly, -10, 10);
        rootsVec.push_back(roots);
        size_t rootNumber = roots.size();
        rootsNumberVec.push_back(rootNumber);
    }

    std::vector<double> allRoots;
    for (const auto& roots : rootsVec) {
        allRoots.insert(allRoots.end(), roots.begin(), roots.end());
    }

    // Calculate and print average number of roots
    double sum = std::accumulate(rootsNumberVec.begin(), rootsNumberVec.end(), 0.0);
    double average = sum / rootsNumberVec.size();
    std::cout << average << std::endl;

    plotresult3a_(allRoots); // Plot roots for exercise 3a
    */
    

}