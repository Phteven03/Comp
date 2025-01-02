#define _USE_MATH_DEFINES
#include <iostream>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>
#include <vectormath.h>

#include "mathfunc3.h"
#include "plot3.h"

const long double M_G = 6.67430e-11;


int main() {
    //-------exercise 1 --------
    /*
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
    */

    //------- exercise 3----------

    std::vector<long double> R1 = { -1.0, 0, 0 };
    std::vector<long double> R2 = { 1.0, 0, 0 };
    //std::vector<long double> r = { 0, std::sqrt(3)/2, 0};
    std::vector<long double> r = { -1.2, 0, 0};
    std::vector<long double> velocity = { 0, 0.001, 0 };
    std::vector<long double> omega = { 0, 0, 0.0001 };
    long double m = 1.0e-3;
    long double M1 = 1.0e5;
    long double M2 = 1.0e2;
    long double dt = 0.0001;
    size_t n = 1000;

    std::vector<long double> totalForce = totalForce_(r, R1, R2, velocity, omega, m, M1, M2);
    std::vector<long double> totalForceDimLess = totalForceDimLess_(r, R1, R2, velocity, omega, m, M1, M2);


    std::vector<std::vector<long double>> ri = forwardEuler_(r, R1, R2, velocity, omega, m, M1, M2, dt, n);
    std::vector<long double> xi;
    std::vector<long double> yi;
    std::vector<long double> zi;
    for (const auto& pos : ri) {
            xi.push_back(pos[0]);
            yi.push_back(pos[1]);
            zi.push_back(pos[2]);
    }
    std::vector<long double> Rx = { R1[0],R2[0] };
    std::vector<long double> Ry = { R1[1],R2[1] };
    std::vector<long double> Rz = { R1[2],R2[2] };

    matplot::plot3(xi,yi,zi);
    matplot::hold(matplot::on);
    matplot::plot3(Rx, Ry, Rz,"r.");
    matplot::hold(matplot::off);
    matplot::xlabel("X");
    matplot::ylabel("Y");
    matplot::zlabel("Z");
    matplot::show();
}