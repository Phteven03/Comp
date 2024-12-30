#include <iostream>
#include <vector>
#include <cmath>


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

	//theta
	for (auto& charge : charges) {

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

			charge.theta += deltaTheta;
			charge.phi += deltaPhi;
		}
	}
}

