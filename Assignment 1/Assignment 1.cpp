#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>
#include "vectormath.h"
#include "matplot/matplot.h"

#ifdef _DEBUG
#define COMMENT_START /*
#define COMMENT_END */
#else 
#define COMMENT_START
#define COMMENT_END
#endif

//exercise 1
#include "mathfunc1a.h"
#include "plot1.h"

//exercise 2
#include "mathfunc2a.h"


const double PI = 3.14159265359;
const double SINGULARITY_ERROR_TERM = 1e-5;

//exercise 1a

static double integrand1a_(double y) {
    return 1 / std::sqrt(1 - std::pow(y, 4));
}
static double trapezoidRule_(double a, double b, int n, std::function<double(double)> integrand) {
    if (isinf(integrand(a))) a += SINGULARITY_ERROR_TERM;
    if (isinf(integrand(b))) b -= SINGULARITY_ERROR_TERM;
    double stepWidth = (b - a) / n;
    double trapIntegralValue = 0.5 * (integrand(a) + integrand(b));
    for (int i = 1; i < n; ++i) {
        trapIntegralValue += integrand(a + i * stepWidth);
    }
    return trapIntegralValue * stepWidth;
}
static double simpsonRule_(double a, double b, int n, std::function<double(double)> integrand) {
    double simpIntegralValue = 0;
    double stepWidth = (b - a) / n;
    double oddSumValue = 0;
    double evenSumValue = 0;

    for (int i = 1; i < n; i += 2) {
        oddSumValue += integrand(a + i * stepWidth);
    }


    for (int i = 2; i <= n; i += 2) {
        evenSumValue += integrand(a + i * stepWidth);
    }

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
    return 1 / std::sqrt((std::exp(a) - std::exp(x)));
}
static double integrand1b3_(double x, double a) {
    return 1 / std::sqrt(std::cos(x) - std::cos(a));
}

//exercise 1c

static double potential_(double k, double x, double d) {
    return std::tanh(k * x);
}
static double integrand1c1_(double x, double k) {
    return 1 / std::sqrt(std::tanh(k) - std::tanh(k * x));
}
static double integrand1c2_(double x, double k) {
    return 1 / std::sqrt(std::tanh(k / 2.0) - std::tanh(k * x));
}


//exercise 2a
double mu = 3e-6;
std::vector<double> function2a = { -mu, 2 * mu, -mu, 3 - 2 * mu, mu - 3, 1 };


int main() {

    //exercise 1a
    /*const double a = 0;
    const double b = 1;
    const double n = 100;


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

    for (int i = 1; i < n; ++i) {
        gaussIntValue = gaussianQudrature_(a, b, i, integrand1a_);
        gaussIntVector.push_back(gaussIntValue);

        trapIntValue = trapezoidRule_(a, b, i, integrand1a_);
        trapIntVector.push_back(trapIntValue);
    }
    for (int i = 1; i < 2 * n; i += 2)
    {
        simpIntValue = simpsonRule_(a , b , i , integrand1a_);
        simpIntVector.push_back(simpIntValue);
    }
    plotresult1a_(trapIntVector, simpIntVector, gaussIntVector, n);*/

    //exercise 1b
    /*std::vector<std::vector<double>> allIntegrals;
    int n1b = 50; 
    for (double a = 0; a <= PI; a += (PI / 1e2)) {

        std::vector<std::function<double(double)>> integrandFixedA = {
            [a](double x) { return integrand1b1_(x, a); },
            [a](double x) { return integrand1b2_(x, a); },
            [a](double x) { return integrand1b3_(x, a); }
        };

        std::vector<double> intVectorForA;

        for (const auto& integrand : integrandFixedA) {
            double intVal = gaussianQudrature_(0, a, n1b , integrand);
            intVectorForA.push_back(intVal);
        }

        allIntegrals.push_back(intVectorForA);
    }

    plotresult1b_(allIntegrals);

    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b1_: " << allIntegrals[i][0] << std::endl;
    }
    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b2_: " << allIntegrals[i][1] << std::endl;
    }
    for (size_t i = 0; i < allIntegrals.size(); ++i) {
        std::cout << "Result for a[" << i << "] with integrand1b3_: " << allIntegrals[i][2] << std::endl;
    }*/

    ////exercise 1c 
    /*std::vector<double> potential1c1Vector;
    std::vector<double> potential1c05Vector;
    for (double k = 0; k <= 20; ++k) {

        auto integrandFixedK = [k](double x) { return integrand1c1_(x, k); };
        double potential1c1 = gaussianQudrature_(0, 1, 50, integrandFixedK);
        potential1c1Vector.push_back(potential1c1);

        auto integrandFixedK05 = [k](double x) { return integrand1c2_(x, k); };
        double potential1c05 = gaussianQudrature_(0, 0.5, 50, integrandFixedK05);
        potential1c05Vector.push_back(potential1c05);
    }
    printVector(potential1c1Vector);
    printVector(potential1c05Vector);*/
    

    //exercise 2a
    std::vector<double> functiontest = { 6, -5, 1 };
    std::vector<double> convergenceNewton = newtonConvergence_(function2a, -3, 3);
    std::vector<double> convergenceBisection = bisectonConvergence_(function2a, -3, 3);

    std::vector<double> rootsbis = polyRootBisection_(function2a, -3, 3);
    printVector(rootsbis);
    std::vector<double> rootsnet = polyRootNewtonRaphson_(function2a, -3, 3);
    printVector(rootsnet);

    //plotresult2a_(convergenceBisection, convergenceNewton);

    //exercise 2b 

    std::vector<double> errorIVec;
    std::vector<double> errorIp1Vec;
    for (size_t i = 0; i < convergenceNewton.size() - 1; ++i) {
        double errorI = std::log(std::abs(convergenceNewton[i] - convergenceNewton[convergenceNewton.size()]));
        errorIVec.push_back(errorI);
        double errorIp1 = std::log(std::abs(convergenceNewton[i + 1] - convergenceNewton[convergenceNewton.size()]));
        errorIp1Vec.push_back(errorIp1);
    }

    //double k = errorIVec[]

    //plotresult2b_(errorIVec, errorIp1Vec);
    



}