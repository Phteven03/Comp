#define _USE_MATH_DEFINES
#include <iostream>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>

#include "mathfunc3.h"
#include "plot3.h"

int main() {
    //-------exercise 1 --------
    const int numCharges = 5;
    const double stepSize = 1e-3;
    const double maxIterations = 100000;


    std::vector<double> chargeX;
    std::vector<double> chargeY;
    std::vector<double> chargeZ;

    std::vector<charge> charges(numCharges);
    for (auto& charge : charges) {
        charge.theta = M_PI * (std::rand() / (double)RAND_MAX);
        charge.phi = 2 * M_PI * (std::rand() / (double)RAND_MAX);
    }

    gradientDecent_(charges, stepSize, maxIterations);

    for (auto& charge : charges) {
        chargeX.push_back(std::sin(charge.theta) * std::cos(charge.phi));
        chargeY.push_back(std::sin(charge.theta) * std::sin(charge.phi));
        chargeZ.push_back(std::cos(charge.theta));
    }

    plotresult1c_(chargeX, chargeY, chargeZ);

}