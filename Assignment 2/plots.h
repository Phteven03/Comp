#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "matplot/matplot.h"
#include "mathfunc2.h"
#include <string>

#include <set>
#include <thread>

void plotResult1b(std::vector<float>& fftTimes, std::vector<float>& dftTimes);
void plotResult1c(std::vector<double>& frequencies, std::vector<double>& powers);
void plotResult2b(splineValues splineValues);
void plotResult2c(std::vector<double>& omega, std::vector<size_t>& iterations);
void plotResult3c(std::vector<double>& z, std::vector<std::vector<double>>& eigenVectorMatrix, std::vector<double>& eigenValues);