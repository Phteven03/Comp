#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector> 

#include "matplot/matplot.h"
#include "plots.h"
#include "mathfunc2.h"


void plotResult1b(std::vector<float>& fftTimes, std::vector<float>& dftTimes) {
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

void plotResult1c(std::vector<double>& frequencies, std::vector<double>& powers) {
    matplot::plot(frequencies, powers, ".");
    matplot::xlabel("f / Hz");
    matplot::ylabel("P / Watt / Hz");
    matplot::legend({ "Powerspectrum" , "0" });
    matplot::title("Power spectral density");
    matplot::grid(matplot::on);
    matplot::xlim({ 0, 500 });
    matplot::show();

}

void plotResult2b(splineValues splineValues) {
    matplot::plot(splineValues.xValues, splineValues.splineValues);
    matplot::xlabel("x");
    matplot::ylabel("y");
    matplot::title("Spline Fit");
    matplot::grid(matplot::on);
    matplot::show();
}

void plotResult2c(std::vector<double>& omega, std::vector<size_t>& iterations) {

    std::vector<int> intIterations;
    for (size_t i = 0; i < iterations.size(); ++i) {
        intIterations.push_back(static_cast<int>(iterations[i]));
    }
    matplot::plot(omega, intIterations);
    matplot::xlabel("Omega");
    matplot::ylabel("Iterations");
    matplot::grid(matplot::on);
    matplot::show();
}

void plotResult3c(std::vector<double>& z, std::vector<std::vector<double>>& eigenVectorMatrix, std::vector<double>& eigenValues) {
    size_t n = z.size();

    std::vector<double> zFirstStrand(z.begin(), z.begin() + n / 2);
    std::vector<double> zSecondStrand(z.begin() + n / 2, z.end());

    std::vector<double> firstStrandEigenValues(eigenValues.begin(), eigenValues.begin() + n / 2);
    std::vector<double> secondStrandEigenValues(eigenValues.begin() + n / 2, eigenValues.end());

    std::vector<std::vector<double>> firstStrandMatrix;
    std::vector<std::vector<double>> secondStrandMatrix;

    for (const auto& row : eigenVectorMatrix) {
        size_t mid = row.size() / 2;
        std::vector<double> firstHalf(row.begin(), row.begin() + mid);
        std::vector<double> secondHalf(row.begin() + mid, row.end());

        firstStrandMatrix.push_back(firstHalf);
        secondStrandMatrix.push_back(secondHalf);
    }



    for (size_t i = 0; i < 10; ++i) {
        matplot::figure(i);
        matplot::plot(zFirstStrand, firstStrandMatrix[i]);
        matplot::hold(matplot::on);
        matplot::plot(zSecondStrand, secondStrandMatrix[i]);
        matplot::hold(matplot::off);
        matplot::xlabel("z");
        matplot::ylabel("eigenvectors");
        matplot::legend({ "strand 1", "strand 2" });
        matplot::grid(matplot::on);
        matplot::show();

    }



}