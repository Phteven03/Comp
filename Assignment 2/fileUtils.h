#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "mathfunc2.h"
#include "matplot/matplot.h"
#include "fftw3.h"

std::vector<std::vector<double>> readTxt2Matrix_(const std::string& fileName);