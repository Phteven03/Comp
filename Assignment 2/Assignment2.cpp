#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <numeric>

#include "mathfunc2.h"
#include "fileUtils.h"
#include "matplot/matplot.h"
#include "fftw3.h"
#include "plots.h"



std::vector<std::complex<double>> discreteFourierTransform_(std::vector<double>& values, StepTimer* stepTimer) {
    if (stepTimer) {
        stepTimer->startTimer();
    }
    size_t n = values.size();
    double D = 2.0 * M_PI / n;
    std::vector<std::complex<double>> fourierTransformed;
    fourierTransformed.reserve(n);

    for (size_t k = 0; k < n; ++k) {
        std::complex<double> sum(0, 0);
        for (int j = 0; j < n; ++j) {
            double realPart = std::cos(k * j * D);
            double imagPart = std::sin(k * j * D);
            std::complex<double> w(realPart, -imagPart);
            sum += w * values[j];
        }
        fourierTransformed.push_back(sum / static_cast<double>(n));
    }
    if (stepTimer) {
        stepTimer->stopStoreTimer();
    }

    return fourierTransformed;
}

int main() {


    //-------- exercise 1 --------------
    std::vector<std::vector<double>> data = readTxt2Matrix_("single_tone.txt");
    std::vector<double> dataLeft = data[0];
    std::vector<double> dataRight = data[1];
    size_t dataSize = dataLeft.size();
    int sampleRate = 44100;
    /*
    //--------exercise 1a ----------------

    
    std::vector<std::complex<double>> fftResult = FFT_(dataLeft);

    std::vector<std::complex<double>> dftResult = discreteFourierTransform_(dataLeft);
    size_t sizefft = fftResult.size();

    std::vector<double> sum(sizefft);

    for (size_t i = 0; i < sizefft; ++i) {
        sum[i] = std::abs(std::abs(fftResult[i]) / (sizefft / 2.0) - std::abs(dftResult[i]) / (sizefft / 2.0));
        if (std::abs(sum[i]) < 1e-7) {
            sum[i] = 0;
        }
    }
    double totalSum = std::accumulate(sum.begin(), sum.end(), 0.0);

    std::cout << "Difference between FFT and DFT: " << totalSum << std::endl;

    //------------- exercise 1b --------------

    StepTimer stepTimerfft;
    StepTimer stepTimerdft;
    std::vector<std::complex<double>> fftResultTimed;
    std::vector<std::complex<double>> dftResultTimed;
    std::vector<double> m;
    for (size_t m = 1; m < std::floor(sizefft/1000); ++m) {
        std::vector<double> dataSubset(dataLeft.begin(), dataLeft.begin() + m);
        fftResultTimed = FFT_(dataSubset, &stepTimerfft);
        dftResultTimed = discreteFourierTransform_(dataSubset, &stepTimerdft);
    }
    std::vector<float> fftTimes = stepTimerfft.getTimes();
    std::vector<float> dftTimes = stepTimerdft.getTimes();

    plotResult1b(fftTimes, dftTimes);
    */

    //-------- exercise 1c --------------

    std::vector<std::pair<double, double>> powerSpectrumData = powerSpectrum_(dataLeft, sampleRate);
    std::vector<double> frequencies;
    std::vector<double> powers;
    for (const auto& pair : powerSpectrumData) {
        frequencies.push_back(pair.first);
        powers.push_back(pair.second);
    }

    matplot::plot(frequencies, powers);
    matplot::show();

}