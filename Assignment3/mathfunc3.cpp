#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <vectormath.h>
#include <functional>
#include "mathfunc2a.h"
#include "mathfunc1a.h"

#include "mathfunc3.h"


//----- Function for exercise 1 

double potential_(double x, double m, double lambda, double omega) {
    return 0.5 * m * m * omega * omega * x * x + (lambda / 24.0) * x * x * x * x;
}

double orthodrome_(const charge& q1, const charge& q2) {
	return std::acos(std::sin(q1.theta) * std::sin(q2.theta) * std::cos(q1.phi - q2.phi) + std::cos(q1.theta) * std::cos(q2.theta));
}

double rij_(const charge& q1, const charge& q2) {
	double x1 = std::sin(q1.theta) * std::cos(q1.phi);
	double y1 = std::sin(q1.theta) * std::sin(q1.phi);
	double z1 = std::cos(q1.theta);

	double x2 = std::sin(q2.theta) * std::cos(q2.phi);
	double y2 = std::sin(q2.theta) * std::sin(q2.phi);
	double z2 = std::cos(q2.theta);

	return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

double totalEnergy_(std::vector<charge> charges) {
	double energy = 0.0;
	size_t vecSize = charges.size();
	for (size_t i = 0; i < vecSize; ++i) {
		for (size_t j = i + 1; j < vecSize; ++j) {
			energy += 1 / rij_(charges[i], charges[j]); // Energy contribution from each pair of charges
		}
	}
	return energy;
}

void numericalGradient_(std::vector<charge>& charges, double stepSize) {
	double currentEnergy = totalEnergy_(charges);

	for (auto& charge : charges) {
		// Compute gradient with respect to theta
		double originalTheta = charge.theta;
		charge.theta += stepSize;
		double energyTheta = totalEnergy_(charges);
		charge.gradTheta = (energyTheta - currentEnergy) / stepSize;
		charge.theta = originalTheta;

		// Compute gradient with respect to phi
		double originalPhi = charge.phi;
		charge.phi += stepSize;
		double energyPhi = totalEnergy_(charges);
		charge.gradPhi = (energyPhi - currentEnergy) / stepSize;
		charge.phi = originalPhi;
	}
}

void gradientDecent_(std::vector<charge>& charges, double stepSize, double maxIterations) {
	for (size_t iter = 0; iter < maxIterations; ++iter) {
		numericalGradient_(charges, 1e-4); // Compute numerical gradients

		for (auto& charge : charges) {
			// Update charge positions using gradient descent
			charge.theta += -stepSize * charge.gradTheta;
			charge.phi += -stepSize * charge.gradPhi;
		}
	}
}


// -------- Functions for exercise 2 

//------ Functions for exercise 3

std::vector<double> gravitationalForce_(std::vector<double>& r, std::vector<double>& R, double massLight, double massHeavy) {
	const double M_G = 6.67430e-11; // Gravitational constant
	std::vector<double> diffVec = r - R;
	double distance = norm_(diffVec);
	return diffVec * (-M_G * massLight * massHeavy / (distance * distance * distance));
}

std::vector<double> inertialForce_(std::vector<double>& r, std::vector<double>& velocity, std::vector<double>& omega, double massLight) {
	// Centrifugal force
	std::vector<double> centrifugalForce = euklidCrossProduct_(omega, euklidCrossProduct_(omega, r)) * massLight;
	// Coriolis force
	std::vector<double> coriolisForce = euklidCrossProduct_(omega, velocity) * (-2.0 * massLight);
	return centrifugalForce + coriolisForce;
}

std::vector<double> totalForce_(std::vector<double>& r, std::vector<double>& R1, std::vector<double>& R2, std::vector<double>& velocity, std::vector<double>& omega, double m, double M1, double M2) {
	// Gravitational forces from both heavy bodies
	std::vector<double> grav1 = gravitationalForce_(r, R1, m, M1);
	std::vector<double> grav2 = gravitationalForce_(r, R2, m, M2);
	// Inertial forces
	std::vector<double> inert = inertialForce_(r, velocity, omega, m);
	return grav1 + grav2 + inert;
}

std::vector<double> totalForceDimLess_(std::vector<double> rNew, std::vector<double> R1New, std::vector<double> R2New, std::vector<double> omegaNew, std::vector<double> vNew, double mu) {
	double dist1 = norm_(rNew - R1New);
	double dist2 = norm_(rNew - R2New);

	// Gravitational forces in the dimensionless system
	std::vector<double> F1 = -mu / (dist1 * dist1 * dist1) * (rNew - R1New);
	std::vector<double> F2 = -(1 - mu) / (dist2 * dist2 * dist2) * (rNew - R2New);
    //std::cout << "Grav:" << std::endl;
    //printVector(F1 + F2);

	// Inertial force
	std::vector<double> FI = -2 * euklidCrossProduct_(omegaNew, vNew) - euklidCrossProduct_(omegaNew, euklidCrossProduct_(omegaNew, rNew));
    //std::cout << "Intertial:" << std::endl;
    //printVector(FI);
	return F1 + F2 + FI;
}

std::vector<std::vector<double>> forwardEuler_(std::vector<double> rNew, std::vector<double> R1New, std::vector<double> R2New, std::vector<double> omegaNew, std::vector<double> vNew, double mu, double dt, int n) {
	std::vector<std::vector<double>> position;
	position.push_back(rNew);

	for (size_t i = 0; i < n; ++i) {
		std::vector<double> a = totalForceDimLess_(rNew, R1New, R2New, omegaNew, vNew, mu);
		vNew = vNew + dt * a;  // Update velocity
		rNew = rNew + dt * vNew;  // Update position
		position.push_back(rNew);
	}
	return position;
}

std::vector<std::vector<double>> lagrangePointFinder_(double mu) {
    std::vector<std::vector<double>> lagrangePoints;

    double muNew = static_cast<double>(mu);
    double leftLimit = -2; // Left limit for root-finding.
    double rightLimit = 2; // Right limit for root-finding.

    // Polynomial coefficients for L1, L2, and L3 points.
    std::vector<std::vector<double>> LPolys = {
        { -muNew, 2 * muNew, -muNew, 3 - 2 * muNew, muNew - 3, 1 },
        { -muNew, -2 * muNew, -muNew, 3 - 2 * muNew, 3 - muNew, 1 },
        { -muNew, 12 + 14 * muNew, -24 - 13 * muNew, 6 * muNew + 19, -7 - muNew, 1 }
    };

    for (const auto& Li : LPolys) {
        // Find roots of the polynomial using Newton-Raphson and adjust for offset
        std::vector<double> roots = polyRootNewtonRaphson_(Li, leftLimit, rightLimit);
        roots.push_back(0.0); // Include z-coordinate for 2D.
        lagrangePoints.push_back(roots);
    }

    // Adjust coordinates for physical positions.
    lagrangePoints[0][0] = (1 - mu) - lagrangePoints[0][0]; // L1 x-coordinate.
    lagrangePoints[1][0] += 1 - mu; // L2 x-coordinate.
    lagrangePoints[2][0] = -1 - lagrangePoints[2][0]; // L3 x-coordinate.

    // Add L4 and L5 coordinates based on trigonometric positions.
    std::vector<double> L4 = { std::cos(M_PI / 3) - muNew, std::sin(M_PI / 3) };
    std::vector<double> L5 = { std::cos(-M_PI / 3) - muNew, std::sin(-M_PI / 3) };
    lagrangePoints.push_back(L4);
    lagrangePoints.push_back(L5);

    return lagrangePoints;
}


// ------ Functions for exercise 4

std::vector<double> rungeKutta_(double x, std::vector<double> y, double h, std::function<std::vector<double>(double, std::vector<double>)> func) {
    std::vector<double> k1 = func(x, y);
    std::vector<double> k2 = func(x + h / 2, y + h / 2 * k1); // Midpoint evaluation.
    std::vector<double> k3 = func(x + h / 2, y + h / 2 * k2); // Midpoint evaluation based on k2.
    std::vector<double> k4 = func(x + h, y + h * k3); // Endpoint evaluation.

    // Return the weighted sum of increments.
    return (y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4));
}

std::vector<double> schroedinger_(double x, double m, double lambda, double omega, double hbar, double E, std::vector<double>& y) {
    std::vector<double> dydx(2);
    dydx[0] = y[1]; // psi' = psi[1].
    dydx[1] = 2 * m / (hbar * hbar) * (potential_(x, m, lambda, omega) - E) * y[0]; // psi''.
    return dydx;
}

double psiL_(double E, double m, double lambda, double omega, double hbar, double L, double h, bool even) {
    std::vector<double> y = even ? std::vector<double>{1.0, 0.0} : std::vector<double>{ 0.0, 1.0 }; // Initial conditions.
    double x = 0.0;

    auto schroedinger = [&m, &lambda, &omega, &hbar, &E](double x, std::vector<double> y) {
        return schroedinger_(x, m, lambda, omega, hbar, E, y);
        };

    // Integrate until x reaches L.
    while (x < L) {
        y = rungeKutta_(x, y, h, schroedinger);
        x += h;
    }

    return y[0]; // Return the wave function at L.
}

std::vector<std::vector<double>> findBracket_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even) {
    std::vector<std::vector<double>> significantIntervals;

    double EStart = 0.4; // Initial energy guess.
    const double stepWidth = 1e-1; // Step size for energy increment.
    double previousPsi = psiL_(EStart, m, lambda, omega, hbar, L, h, even);
    double previousE = EStart;

    size_t iteration = 0;

    while (iteration < (N / 2)) {
        double currentE = previousE + stepWidth;
        double currentPsi = psiL_(currentE, m, lambda, omega, hbar, L, h, even);

        // Check for a sign change in the wave function.
        if ((previousPsi * currentPsi) < 0.0) {
            significantIntervals.push_back({ previousE, currentE });
            iteration++;
        }
        previousPsi = currentPsi;
        previousE = currentE;
    }

    return significantIntervals;
}

std::vector<double> eigenEnergyFinder_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even, double tolerance) {
    // Find initial energy brackets where roots exist.
    std::vector<std::vector<double>> guessIntervals = findBracket_(m, lambda, omega, hbar, L, h, N, even);

    std::vector<double> energies;

    for (const auto& interval : guessIntervals) {
        double leftLimit = interval[0];
        double rightLimit = interval[1];
        double leftValue = psiL_(leftLimit, m, lambda, omega, hbar, L, h, even);
        double rightValue = psiL_(rightLimit, m, lambda, omega, hbar, L, h, even);
        double midpoint;

        // Apply bisection method to find eigenenergies.
        while ((rightLimit - leftLimit) / 2.0 > tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            double midpointValue = psiL_(midpoint, m, lambda, omega, hbar, L, h, even);

            if (leftValue * midpointValue < 0) { // Root is in the left subinterval.
                rightLimit = midpoint;
                rightValue = midpointValue;
            }
            else { // Root is in the right subinterval.
                leftLimit = midpoint;
                leftValue = midpointValue;
            }
        }
        energies.push_back(midpoint); // Add the converged eigenenergy.
    }

    return energies;
}


