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
    const int numCharges = 6;
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
    
    //start
    std::vector<long double> R1 = { -1.0, 0, 0 };
    std::vector<long double> R2 = { 1.0, 0, 0 };
    //std::vector<long double> r = { 0, std::sqrt(3)/2, 0};
    std::vector<long double> r = { 0, 1, 0};
    std::vector<long double> velocity = { 0, -0.1, 0 };
    std::vector<long double> omega = { 0, 0, -1000 };
    long double M1 = 1.0e3;
    long double M2 = 1.0e5;
    long double dt = 0.001;
    size_t n = 10000;

    //scaling factors
    long double L = norm_(R2 - R1);
    long double M = M1 + M2;
    long double T = 1 / std::abs(omega[2]);

    //funciton variabls dimensionless
    std::vector<long double> rNew = 1 / L * r;
    std::vector<long double> R1New = 1 / L * R1;
    std::vector<long double> R2New = 1 / L * R2;
    std::vector<long double> vNew = (T / L) * velocity;
    std::vector<long double> omegaNew = omega * T;
    long double mu = M2 / M;

    std::vector<std::vector<long double>> ri = forwardEuler_(rNew,R1New,R2New,omegaNew,vNew,mu,dt,n);



    std::vector<long double> Rx = { -(1-mu),mu};
    std::vector<long double> Ry = { 0,0 };


    std::vector<long double> xi;
    std::vector<long double> yi;
    std::vector<long double> zi;
    for (const auto& pos : ri) {
            xi.push_back(pos[0]);
            yi.push_back(pos[1]);
            zi.push_back(pos[2]);
    }

    matplot::plot(xi,yi);
    matplot::hold(matplot::on);
    matplot::plot(Rx, Ry, "r.");
    matplot::text(Rx[0], Ry[0], "M1");
    matplot::text(Rx[1], Ry[1], "M2");
    matplot::xlabel("X");
    matplot::ylabel("Y");
    matplot::zlabel("Z");
    matplot::hold(matplot::off);
    matplot::show(); 
    
}