#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "matplot/matplot.h"
#include "plots.h"

void plotResult1b(std::vector<float> fftTimes, std::vector<float> dftTimes) {
    matplot::plot(fftTimes);
    matplot::hold(matplot::on);
    matplot::plot(dftTimes);
    matplot::xlabel("m");
    matplot::ylabel("t / ms");
    matplot::legend({ "FFT-Time", "DFT-Time" });
    matplot::title("Computation Time FFT vs. DFT");
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);
    matplot::show();

}

void plotResult1c(std::vector<double> frequencies, std::vector<double> powers) {
    matplot::plot(frequencies, powers);
    matplot::xlabel("f / Hz");
    matplot::ylabel("P / Watt / Hz");
    matplot::legend({ "Powerspectrum" , "0"});
    matplot::title("Power spectral density");
    matplot::grid(matplot::on);
    matplot::xlim({ 0, 500 });
    matplot::show();

}