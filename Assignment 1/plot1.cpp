#include <iostream>
#include <vector>
#include "plot1.h"
#include "vectormath.h"
#include "matplot/matplot.h"
#include "mathfunc1a.h"

void plotresult1a_(std::vector<double> trapIntVector, std::vector<double> simpIntVector, std::vector<double> gaussIntVector,  int n) {
	const double REALINTVALUE = 1.311028777146120;
	std::vector<double> nVector;
	for (int j = 0; j <= n; ++j) {
		nVector.push_back(j);
	}

	size_t sizeIntVector = trapIntVector.size();
	for (size_t i = 0; i < sizeIntVector; ++i) {
		trapIntVector[i] = abs(trapIntVector[i] - REALINTVALUE);
		gaussIntVector[i] = abs(gaussIntVector[i] - REALINTVALUE);
	}
	
	size_t simpSizeIntVector = simpIntVector.size();
	for (size_t i = 0; i < simpSizeIntVector; ++i) {
		simpIntVector[i] = abs(simpIntVector[i] - REALINTVALUE);
	}

	matplot::semilogy(nVector, gaussIntVector);
	matplot::hold(matplot::on);
	matplot::semilogy(nVector, trapIntVector);
	matplot::semilogy(nVector, simpIntVector);
	matplot::xlabel("N");
	matplot::ylabel("Error");
	matplot::legend({ "Gaussian Quadrature", "Trapazoid", "Simpson" });
	matplot::grid(matplot::off);
	matplot::title("Error of different integration methods");
	matplot::show();
}

void plotresult1b_(std::vector<std::vector<double>> allIntegrals) {
	size_t size = allIntegrals.size();
	std::vector<double> nVector(size);
	std::vector<double> yVector(size);
	for (int j = 0; j < 3; ++j) {
		for (size_t i = 0; i < size; ++i) {
			nVector[i] = i;
			yVector[i] = allIntegrals[i][j];
		}
		matplot::plot(nVector, yVector);
		matplot::hold(matplot::on);
	}

	matplot::xlabel("a");
	matplot::ylabel("T(a)");
	auto legend = matplot::legend({ "V(x) = cosh(x)", "V(x) = exp(|x|)", "V(x) = -cos(x)" }); 
	legend->location(matplot::legend::general_alignment::topleft);
	matplot::title("differenz V(x)");
	matplot::hold(matplot::off);
	matplot::show();
}

void plotresult2a_(std::vector<double> convergenceBisection, std::vector<double> convergenceNewton) {
	size_t sizeBi = convergenceBisection.size();
	std::vector<double> nVectorBi;
	for (size_t i = 0; i < sizeBi; ++i) {
		nVectorBi.push_back(i);
	}

	size_t sizeNe = convergenceNewton.size();
	std::vector<double> nVectorNe;
	for (size_t j = 0; j < sizeNe; ++j) {
		nVectorNe.push_back(j);
	}

	matplot::plot(nVectorBi, convergenceBisection);
	matplot::hold(matplot::on);
	matplot::plot(nVectorNe, convergenceNewton);
	matplot::hold(matplot::off);
	matplot::legend({ "Bisection" , "Newton" });
	matplot::ylabel("root Value");
	matplot::xlabel("Iteration");
	matplot::title("Convergence");
	matplot::grid(matplot::on);
	matplot::show();
}

void plotresult2b_(std::vector<double> errorI, std::vector<double> errorIp1) {

	matplot::loglog(errorIp1, errorI);
	matplot::grid(matplot::on);
	matplot::xlabel("error i+1");
	matplot::ylabel("error i");
	matplot::show();
}

void plotresult3a_(std::vector<double> allRoots) {

	matplot::hist(allRoots, 500);
	matplot::title("Distribution of roots");
	matplot::xlabel("Location");
	matplot::ylabel("Number of Roots");
	matplot::grid(matplot::on);
	matplot::show();
}