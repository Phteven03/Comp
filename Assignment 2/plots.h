#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "matplot/matplot.h"
#include "mathfunc2.h"

void plotResult1b(std::vector<float>& fftTimes, std::vector<float>& dftTimes);
void plotResult1c(std::vector<double>& frequencies, std::vector<double>& powers);
void plotResult2b(splineValues splineValues);