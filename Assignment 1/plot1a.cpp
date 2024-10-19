#include <iostream>
#include <vector>
#include "plot1a.h"
#include "vectormath.h"
#include "matplot/matplot.h"
#include "mathfunc1a.h"

void plotresult1a_(std::vector<double> trapIntVector, std::vector<double> simpIntVector, std::vector<double> gaussIntVector,  int n) {
	std::vector<double> nVector;

	for (int j = 0; j < n; ++j) {
		nVector.push_back(j);
	}

	//matplot::plot(nVector, gaussIntVector);
	//matplot::hold(matplot::on);
	//matplot::plot(nVector, trapIntVector);
	//matplot::plot(nVector, simpIntVector);
	//matplot::xlabel("N");
	//matplot::ylabel("Error");
	//matplot::show();
}

