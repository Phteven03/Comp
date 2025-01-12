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
    //-------exercise 1 & 2 --------
 /*
    const int numCharges = 12;
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

    for (auto& charge : charges) {
        chargeX.push_back(std::sin(charge.theta) * std::cos(charge.phi));
        chargeY.push_back(std::sin(charge.theta) * std::sin(charge.phi));
        chargeZ.push_back(std::cos(charge.theta));
    }

    plotresult1c_(chargeX, chargeY, chargeZ);
    
    // ------- exercise 2 ------ 
    std::vector<charge> charges2(numCharges);
    for (auto& charge2 : charges2) {
        charge2.theta = distTheta(gen);
        charge2.phi = distPhi(gen);
    }
  */  
    
    
    //------- exercise 3----------
    /*
    //L1 = 0.848624 0
    //L2 = 1.14632 0
    //L3 = -1.00082 0
    //L4 = 0.490099 0.866025
    //L5 = 0.490099 - 0.866025
    
    //start
    std::vector<double> R1 = { -1, 0, 0 };
    std::vector<double> R2 = { 1, 0, 0 };
    //std::vector<double> r = { 0, std::sqrt(3)/2, 0};
    std::vector<double> r = { -1, 0.5, 0.0};
    std::vector<double> velocity = { 0, 0, 0 };
    std::vector<double> omega = { 0, 0, 420 };
    double M1 = 1.0e7;
    double M2 = 1.0e5;
    double dt = 0.0001;
    size_t n = 20000;

    //scaling factors
    double M = M1 + M2;
    double T = 2 * M_PI / std::abs(omega[2]);
    double L = norm_(R1 - R2);



    //funciton variabls dimensionless
    double mu = M2 / M;
    std::vector<double> rNew = r;
    std::vector<double> R1New = {-mu,0,0};
    std::vector<double> R2New = {1 - mu,0,0};
    std::vector<double> vNew = velocity * (1 / std::sqrt(M1 + M2));
    printVector(vNew);
    std::vector<double> omegaNew = omega * (1 / std::sqrt(M1 + M2));


    std::vector<std::vector<double>> ri = forwardEuler_(rNew,R1New,R2New,omegaNew,vNew,mu,dt,n);
    std::vector<std::vector<double>> lagrangePoints = lagrangePointFinder_(mu);

    std::cout << "Lagrange Points: " << std::endl;
    printMatrix(lagrangePoints);

    plotresult3b_(ri, lagrangePoints, mu);
    */


    //----- exercise 3 new ----------
    
    std::vector<double> R1 = { -1, 0, 0 };
    std::vector<double> R2 = { 1, 0, 0 };
    std::vector<double> r = { 0.0, 1.0, 0.0 };

    std::vector<double> velocity = { 0, 0, 0 };
    std::vector<double> omega = { 0, 0, 0.01 };

    double M1 = 1.0e22;
    double M2 = 1.0e22;
    double m = 1.0e3;

    std::vector<double> totalforce = totalForce_(r, R1, R2, velocity, omega, m, M1, M2);


    







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