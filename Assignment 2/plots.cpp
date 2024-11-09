#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "matplot/matplot.h"
#include "plots.h"

void plotResult1b(std::vector<float> fftTimes, std::vector<float> dftTimes) {
    matplot::plot(fftTimes);
    matplot::hold(matplot::on);
    matplot::plot(dftTimes);
    matplot::hold(matplot::off);
    matplot::xlabel("m");
    matplot::ylabel("t / ms");
    matplot::legend({ "FFT-Time", "DFT-Time" });
    matplot::title("Computation Time FFT vs. DFT");
    matplot::show();

}