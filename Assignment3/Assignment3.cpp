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

const long double M_G = 6.67430e-11;


static double potential_(double x, double m, double lambda, double omega) {
    return 0.5 * m * m * omega * omega * x * x + (lambda / 24.0) * x * x * x * x;
}

static std::vector<double> schroedinger_(double x, double m, double lambda, double omega, double hbar,  double E, std::vector<double>& y) {
    std::vector<double> dydx(2);
    dydx[0] = y[1];
    dydx[1] = 2 * m / (hbar * hbar) * (potential_(x, m, lambda, omega) - E) * y[0];
    return dydx;
}

static double psiL_(double E, double m, double lambda, double omega, double hbar, double L, double h, bool even) {

    std::vector<double> y = even ? std::vector<double>{1.0, 0.0} : std::vector<double>{ 0.0, 1.0 };
    double x = 0.0;

    auto schroedinger = [&m, &lambda, &omega, &hbar, &E](double x, std::vector<double> y) {
        return schroedinger_(x, m, lambda, omega, hbar, E, y);
    };

    while (x < L) {
        y = rungeKutta_(x, y, h, schroedinger);
        x += h;
    }

    return y[0];
}

static std::vector<std::vector<double>> findBracket_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even) {
    std::vector<std::vector<double>> significantIntervalls;

    double EStart = 0.4;
    const double stepWidth = 1e-1;
    double previousPsi = psiL_(EStart, m, lambda, omega, hbar, L, h, even);
    double previousE = EStart;

    size_t iteration = 0;

    while (iteration < (N / 2)) {
        double currentE = previousE + stepWidth;
        double currentPsi = psiL_(currentE, m, lambda, omega, hbar, L, h, even);

        if ((previousPsi * currentPsi) < 0.0) {
            significantIntervalls.push_back({ previousE, currentE });
            iteration++;
        }
        previousPsi = currentPsi;
        previousE = currentE;
    }

    return significantIntervalls;
}

static std::vector<double> eigenEnergyFinder_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even, double tolerance) {

    // Bracketing
    std::vector<std::vector<double>> guessIntervalls = findBracket_(m,lambda,omega,hbar,L,h,N,even);

    // Bisection-Methode
    std::vector<double> energies;

    for (const auto& interval : guessIntervalls) {
        double leftLimit = interval[0];
        double rightLimit = interval[1];
        double leftValue = psiL_(leftLimit, m, lambda, omega, hbar, L, h, even);
        double rightValue = psiL_(rightLimit, m, lambda, omega, hbar, L, h, even);
        double midpoint;

        // Bisection method to narrow down root intervals
        while ((rightLimit - leftLimit) / 2.0 > tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            double midpointValue = psiL_(midpoint, m, lambda, omega, hbar, L, h, even);

            // Determine which side of the interval contains the root
            if (leftValue * midpointValue < 0) {
                rightLimit = midpoint;
                rightValue = midpointValue;
            }
            else {
                leftLimit = midpoint;
                leftValue = midpointValue;
            }
        }
        energies.push_back(midpoint);
    }

    return energies;
}



int main() {
    //-------exercise 1 --------
    /*
    const int numCharges = 4;
    const double stepSize = 1e-3;
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

    for (auto& charge : charges) {
        chargeX.push_back(std::sin(charge.theta) * std::cos(charge.phi));
        chargeY.push_back(std::sin(charge.theta) * std::sin(charge.phi));
        chargeZ.push_back(std::cos(charge.theta));
    }

    plotresult1c_(chargeX, chargeY, chargeZ);
    */
    
    //------- exercise 3----------
    /*
    //L1 = 0.848624 0
    //L2 = 1.14632 0
    //L3 = -1.00082 0
    //L4 = 0.490099 0.866025
    //L5 = 0.490099 - 0.866025
    
    //start
    std::vector<long double> R1 = { -1, 0, 0 };
    std::vector<long double> R2 = { 1, 0, 0 };
    //std::vector<long double> r = { 0, std::sqrt(3)/2, 0};
    std::vector<long double> r = { 0.848624, 0.0, 0.0};
    std::vector<long double> velocity = { 0, 0, 0 };
    std::vector<long double> omega = { 0, 0, 400 };
    long double M1 = 1.0e7;
    long double M2 = 1.0e5;
    long double dt = 0.0004;
    size_t n = 1000;

    //scaling factors
    //long double L = norm_(R2 - R1);
    long double L = 1;
    long double M = M1 + M2;
    long double T = 1 / std::abs(omega[2]);



    //funciton variabls dimensionless
    long double mu = M2 / M;
    std::vector<long double> rNew = r;
    //std::vector<long double> R1New = 1 / L * R1;
    std::vector<long double> R1New = {-mu,0,0};
    //std::vector<long double> R2New = 1 / L * R2;
    std::vector<long double> R2New = {1 - mu,0,0};
    //std::vector<long double> vNew = (T / L) * velocity;
    std::vector<long double> vNew = velocity * (1 / std::sqrt(M1+M2));
    std::vector<long double> omegaNew = omega * (1 / std::sqrt(M1 + M2));


    std::vector<std::vector<long double>> ri = forwardEuler_(rNew,R1New,R2New,omegaNew,vNew,mu,dt,n);
    std::vector<std::vector<double>> lagrangePoints = lagrangePointFinder_(mu);

    for (const auto& row : lagrangePoints) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

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

    std::cout << "even:" << std::endl;
    //print(evenEigenEnergy.size());
    printVector(evenEigenEnergy);
    std::cout << "odd:" << std::endl;
    //print(oddEigenEnergy.size());
    printVector(oddEigenEnergy);
    */
}