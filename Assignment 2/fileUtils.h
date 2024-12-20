#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>

#include "mathfunc2.h"
#include "matplot/matplot.h"
#include "fftw3.h"

std::vector<std::vector<double>> readTxt2Matrix_(const std::string& fileName);