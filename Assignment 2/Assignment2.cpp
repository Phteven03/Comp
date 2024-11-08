#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>

#include "mathfunc2.h"
#include "fileUtils.h"
#include "matplot/matplot.h"
#include "fftw3.h"



static std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values) {
    Timer time;
    std::vector<std::complex<double>> fourierTransformed;
    int n = values.size();
    double D = 2.0 * M_PI / n;

    for (int k = 0; k < n; ++k) {
        std::complex<double> sum(0, 0);
        for (int j = 0; j < n; ++j) {
            double realPart = std::cos(k * j * D);
            double imagPart = std::sin(k * j * D);
            std::complex<double> w(realPart, -imagPart);
            sum += w * values[j];
        }
        fourierTransformed.push_back(sum / static_cast<double>(n));
    }
    return fourierTransformed;
}

int main() {

    std::vector<std::vector<double>> data = readTxt2Matrix_("single_tone.txt");
    
    size_t n = data.size();
    std::vector<double> xVector;
    std::vector<double> yVector;
    for (size_t i = 0; i < n; ++i) {
        xVector.push_back(i);
    }

    std::vector<std::complex<double>> fftResultlib = FFT_(data[0]);
    std::vector<std::complex<double>> fftResultdisc = discreteFourierTransform_(data[0]);

    //printVector(fftResultdisc);
    //printVector(fftResultlib);
    
    /*
    std::vector<double> fftlibVec;
    std::vector<double> fftdiscVec;

    for (int i = 0; i < n / 2; ++i) {
        double amplitudefftlib = std::abs(fftResultlib[i]) / (n / 2.0);
        fftlibVec.push_back(amplitudefftlib);
        double amplitudefftdisc = std::abs(fftResultdisc[i]) / (n / 2.0);
        fftdiscVec.push_back(amplitudefftdisc);
    }

    matplot::plot(xVector, fftlibVec);
    matplot::hold(matplot::on);
    matplot::plot(xVector, fftdiscVec);
    matplot::hold(matplot::off);
    matplot::show();
    */
}
