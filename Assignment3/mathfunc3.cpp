#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <vectormath.h>
#include "mathfunc2a.h"

#include "mathfunc3.h"



double orthodrome_(const charge& q1, const charge& q2) {
	return std::acos(std::sin(q1.theta) * std::sin(q2.theta) * std::cos(q1.phi - q2.phi) + std::cos(q1.theta) * std::cos(q2.theta));
}

double totalEnergy_(std::vector<charge> charges) {
	double energy = 0.0;
	size_t vecSize = charges.size();
	for (size_t i = 0; i < vecSize; ++i) {
		for (size_t j = i + 1; j < vecSize; ++j) {
			energy += 1 / orthodrome_(charges[i], charges[j]);
		}
	}
	return energy;
}

void numericalGradient_(std::vector<charge>& charges, double stepSize) {
	double currentEnergy = totalEnergy_(charges);

	for (auto& charge : charges) {

		//theta
		double originalTheta = charge.theta;
		charge.theta += stepSize;
		double energyTheta = totalEnergy_(charges);
		charge.gradTheta = (energyTheta - currentEnergy) / stepSize;
		charge.theta = originalTheta;

		//phi
		double originalPhi = charge.phi;
		charge.phi += stepSize;
		double energyPhi = totalEnergy_(charges);
		charge.gradPhi = (energyPhi - currentEnergy) / stepSize;
		charge.phi = originalPhi;
	}
}

void gradientDecent_(std::vector<charge>& charges, double stepSize, double maxIterations) {
	for (size_t iter = 0; iter < maxIterations; ++iter) {

		numericalGradient_(charges, 1e-4);
		for (auto& charge : charges) {
			double deltaTheta = -stepSize * charge.gradTheta;
			double deltaPhi = -stepSize * charge.gradPhi;
			//std::cout << "gradTheta" << std::endl;
			//std::cout << charge.gradTheta << std::endl;
			//std::cout << "gradPhi" << std::endl;
			//std::cout << charge.gradPhi << std::endl;
			charge.theta += deltaTheta; 
			charge.phi += deltaPhi;
		}
	}
}

std::vector<long double> gravitationalForce_(std::vector<long double>& r, std::vector<long double>& R, long double massLight, long double massHeavy) {
	const long double M_G = 6.67430e-11;
	std::vector<long double> diffVec = r - R;
	long double distance = norm_(diffVec);
	return diffVec * (-M_G * massLight * massHeavy / (distance * distance * distance));
}

std::vector<long double> inertialForce_(std::vector<long double>& r, std::vector<long double>& velocity, std::vector<long double>& omega, long double massLight) {
	std::vector<long double> centrifugalForce = euklidCrossProduct_(omega, euklidCrossProduct_(omega, r)) * massLight;
	std::vector<long double> coriolisForce = euklidCrossProduct_(omega, velocity) * (-2.0 * massLight);
	return centrifugalForce + coriolisForce;
}

std::vector<long double> totalForce_(std::vector<long double>& r, std::vector<long double>& R1, std::vector<long double>& R2, std::vector<long double>& velocity, std::vector<long double>& omega, long double m, long double M1, long double M2) {
	std::vector<long double> grav1 = gravitationalForce_(r, R1, m, M1);
	std::vector<long double> grav2 = gravitationalForce_(r, R2, m, M2);
	std::vector<long double> inert = inertialForce_(r, velocity, omega, m);
	return grav1 + grav2 + inert;
}

std::vector<long double> totalForceDimLess_(std::vector<long double> rNew, std::vector<long double> R1New, std::vector<long double> R2New, std::vector<long double> omegaNew, std::vector<long double> vNew, long double mu) {

	long double dist1 = norm_(rNew - R1New);
	long double dist2 = norm_(rNew - R2New);

	//gravitational Forces F1 & F2
	std::vector<long double> F1 = - mu / (dist1 * dist1 * dist1) * (rNew - R1New);
	std::vector<long double> F2 = - (1-mu) / (dist2 * dist2 * dist2) * (rNew - R2New);
	//std::cout << "Grav" << std::endl;
	//printVector(F1 + F2);

	//inertial Force FI
	std::vector<long double> FI = -2 * euklidCrossProduct_(omegaNew, vNew) - euklidCrossProduct_(omegaNew, euklidCrossProduct_(omegaNew, rNew));
	//std::cout << "Inertial" << std::endl;
	//printVector(FI);
	return F1 + F2 + FI;
}

std::vector<std::vector<long double>> forwardEuler_(std::vector<long double> rNew, std::vector<long double> R1New, std::vector<long double> R2New, std::vector<long double> omegaNew, std::vector<long double> vNew, long double mu, long double dt, int n) {
	std::vector<std::vector<long double>> position;

	position.push_back(rNew);
	for (size_t i = 0; i < n; ++i) {
		std::vector<long double> a = totalForceDimLess_(rNew,R1New,R2New,omegaNew,vNew,mu);
		vNew = vNew + dt * a;
		rNew = rNew + dt * vNew;
		position.push_back(rNew);
	}
	return position;
}

std::vector<std::vector<double>> lagrangePointFinder_(long double mu) {
	std::vector<std::vector<double>> lagrangePoints;

	double muNew = static_cast<double>(mu);
	double muNew_ = 1-muNew;
	double leftLimit = -1;
	double rightLimit = 1;

	std::vector<std::vector<double>> LPolys = { { -muNew, 2 * muNew, -muNew, 3 - 2 * muNew, muNew - 3, 1 } , { -muNew, -2 * muNew, -muNew, 3 - 2 * muNew, 3 - muNew, 1 }, { -muNew, 12 + 14 * muNew, -24 - 13 * muNew, 6 * muNew + 19, -7 - muNew, 1 } };

	std::vector<double> muNewVec = { muNew };

	for (const auto& Li : LPolys) {
		std::vector<double> roots = polyRootNewtonRaphson_(Li, leftLimit, rightLimit)-muNewVec;
		roots.push_back(0.0);
		printVector(roots);
		lagrangePoints.push_back(roots);
	}

	double r = 1 + muNew_ * muNew_ - 2 * muNew_ * std::cos(M_PI / 3);
	double phi = std::acos((r * r + muNew_ * muNew_ - 1) / (2 * r * muNew_));
	print(phi);
	std::vector<double> L4 = {r * std::cos(phi) , r * std::sin(phi) };
	std::vector<double> L5 = {r * std::cos(-phi) , r * std::sin(-phi) };

	lagrangePoints.push_back(L4);
	lagrangePoints.push_back(L5);

	return lagrangePoints;
}

