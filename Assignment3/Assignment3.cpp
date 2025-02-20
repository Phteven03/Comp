#define _USE_MATH_DEFINES
#include <iostream>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>
#include <vectormath.h>
#include <random>
#include <functional>
#include "mathfunc2a.h"
#include "mathfunc1a.h"

#include "mathfunc3.h"
#include "plot3.h"

const double M_G = 6.67430e-11;


int main() {
    //-------exercise 1 --------
    StepTimer timer;
	timer.startTimer();
    const int numCharges = 3;
    const double stepSize = 1e-4;
    const double maxIterations = 100000;


    std::vector<double> chargeX;
    std::vector<double> chargeY;
    std::vector<double> chargeZ;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distTheta(0, M_PI);
    std::uniform_real_distribution<> distPhi(0, 2 * M_PI);

    std::vector<charge> charges(numCharges);
    for (auto& charge : charges) {
        charge.theta = distTheta(gen);
        charge.phi = distPhi(gen);
    }

    gradientDecent_(charges, stepSize, maxIterations);
	timer.stopStoreTimer();
	std::vector<float> times = timer.getTimes();
    printVector(times);


    for (auto& charge : charges) {
        chargeX.push_back(std::sin(charge.theta) * std::cos(charge.phi));
        chargeY.push_back(std::sin(charge.theta) * std::sin(charge.phi));
        chargeZ.push_back(std::cos(charge.theta));
    }


    plotresult1c_(chargeX, chargeY, chargeZ);
    
    //----- exercise 3 ----------
    /*
    //0.848624 0
    //1.14632 0
    //-1.00082 0
    //0.490099 0.866025
    //0.490099 - 0.866025



    std::vector<double> r = { 1.0, 1.0, 0.0 };

    std::vector<double> velocity = { 0.0, 0.0, 0.0 };
    std::vector<double> omega = { 0.0, 0.0, 1 };

    double M1 = 1e24;
    double M2 = 1e22;

	double mu = M2 / (M1 + M2);

    std::vector<double> R1 = { -mu, 0, 0 };
    std::vector<double> R2 = { 1 - mu, 0, 0 };

	double dt = 1e-3;
	double n = 20000;
	std::vector<double> totalForce = totalForceDimLess_(r, R1, R2, omega, velocity, mu);
	printVector(totalForce);
    std::vector<std::vector<double>> ri = rungeKutta4_(r, R1, R2, omega, velocity, mu, dt, n);
    std::vector<std::vector<double>> lagrangePoints = lagrangePointFinder_(mu);

    plotresult3b_(ri, lagrangePoints, mu);
    */

    //------exercise 4 ---------
    /*
    const double m = 1.0;
    const double lambda = 24;
    const double omega = 1.0;
    const double hbar = 1.0;
    const double L = 10.0;
    const double h = 0.01; //stepSize for Runge-Kutta
    const double tolerance = 1e-10;
    const double N = 10;

    std::vector<double> evenEigenEnergy = eigenEnergyFinder_(m,  lambda,  omega,  hbar,  L,  h,  N, true ,  tolerance);
    std::vector<double> oddEigenEnergy = eigenEnergyFinder_(m, lambda, omega, hbar, L, h, N, false, tolerance);

    std::cout << "even Eigenenergies:" << std::endl;
    //print(evenEigenEnergy.size());
    printVector(evenEigenEnergy);
    std::cout << "odd Eigenenergies:" << std::endl;
    //print(oddEigenEnergy.size());
    printVector(oddEigenEnergy);
    */
}