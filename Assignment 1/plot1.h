#pragma once
#include <iostream>
#include <vector>
#include "plot1.h"
#include "vectormath.h"
#include "matplot/matplot.h"
#include "mathfunc1a.h"

void plotresult1a_(std::vector<double> trapIntVector, std::vector<double> simpIntVector, std::vector<double> gaussIntVector, int n);

void plotresult1b_(std::vector<std::vector<double>> allIntegrals);

void plotresult2a_(std::vector<double> rootsbisection, std::vector<double> rootsnewton);

void plotresult2b_(std::vector<double> errorI, std::vector<double> errorIp1);