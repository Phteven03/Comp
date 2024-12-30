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
#include <algorithm>
#include <regex>
#include <sstream>
#include <vectormath.h>
#include <iterator>
#include <thread>

#include "mathfunc2.h"
#include "fileUtils.h"
#include "matplot/matplot.h"
#include "fftw3.h"
#include "plots.h"
#include "TridiagonalMatrix.h"
#include "matrixMath.h"



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
        fourierTransformed.push_back(sum);
    }
    if (stepTimer) {
        stepTimer->stopStoreTimer();
    }

    return fourierTransformed;
}

int main() {


    //-------- exercise 1 ----------------
    /*
    std::vector<std::vector<double>> data = readTxt2Matrix_("single_tone.txt");

    std::vector<double> dataLeft = data[0];
    std::vector<double> dataRight = data[1];
    size_t dataSize = dataLeft.size();
    //dataLeft.resize(dataSize + dataSize, 0.0);
    //dataRight.resize(dataSize + dataSize, 0.0);
    //dataSize = dataLeft.size();
    int sampleRate = 44100;
    */

    //-------- exercise 1a ---------------
    /*
    std::vector<std::complex<double>> fftResult = FFT_(dataLeft, -1);
    
    std::vector<std::complex<double>> dftResult = discreteFourierTransform_(dataLeft);


    double sumreal = 0;
    double sumimag = 0;
    for (size_t i = 0; i < dataSize; i++) {
       double diffreal = fftResult[i].real() - dftResult[i].real();
       double diffimag = fftResult[i].imag() - dftResult[i].imag();
       if (std::abs(diffreal) > 1e-12) {
           sumreal += diffreal;
       }
       if (std::abs(diffimag) > 1e-12) {
           sumimag += diffimag;
       }
    }
    std::cout << "Difference between FFT and DFT: " << sumreal << " + i " << sumimag << std::endl;
    */

    //------------- exercise 1b --------------
    /*
    StepTimer stepTimerfft;
    StepTimer stepTimerdft;
    std::vector<std::complex<double>> fftResultTimed;
    std::vector<std::complex<double>> dftResultTimed;

    for (double i = 1; i <= 5e2; ++i) {;
        std::vector<double> partialDataLeft(dataLeft.begin(), dataLeft.begin() + i);
        dftResultTimed = discreteFourierTransform_(partialDataLeft, &stepTimerdft);
        fftResultTimed = FFT_(partialDataLeft, -1, &stepTimerfft);

    }
    std::vector<float> fftTimes = stepTimerfft.getTimes();
    std::vector<float> dftTimes = stepTimerdft.getTimes();

    plotResult1b(fftTimes, dftTimes);
    */

    //-------- exercise 1cde --------------
    /*
    std::vector<std::pair<double, double>> powerSpectrumData = powerSpectrum_(dataLeft, sampleRate);
    std::vector<double> frequencies;
    std::vector<double> powers;
    for (const auto& pair : powerSpectrumData) {
        frequencies.push_back(pair.first);
        powers.push_back(pair.second);
    }

    plotResult1c(frequencies, powers);
    
    //--------- exercise 1de ----------------
    std::vector<double> maxima = maxFinder_(powers);
    std::vector<double> powersSorted = bubbleSort_(maxima);
    size_t n = powersSorted.size();

    std::vector<double> corrPowerValues(dataSize / 2, 0.0);
    std::vector<double> corrFrequencyValue(dataSize / 2);

    for (size_t i = 0; i < powersSorted.size(); ++i) {
        auto it = std::find(powers.begin(), powers.end(), powersSorted[i]);
        size_t index = std::distance(powers.begin(), it);
        corrPowerValues[index] = powers[index];
        corrFrequencyValue[i] = frequencies[index];
    }
    //print(frequencies.size());
    //print(corrPowerValues.size());
    //matplot::figure(1);
    //    matplot::plot(frequencies, corrPowerValues);
    //    matplot::show();
    std::cout << "Fundamental Frequenz: " << corrFrequencyValue[0] << std::endl;
    std::cout << "The note is a D_3" << std::endl;

    std::vector<std::complex<double>> invfftidealValues = FFT_(corrPowerValues, 1);
    size_t sizeInvFFT = invfftidealValues.size();
    std::vector<double> realInvFFT(sizeInvFFT);
    for (size_t i = 0; i < sizeInvFFT; ++i) {
        realInvFFT[i] = std::real(invfftidealValues[i]);
    }

    //matplot::plot(realInvFFT);
    //matplot::show();
    */

    //--------- exercise 2abc --------------
    /*
    std::vector<double> x = { -5,-4,-3,-2,-1,0,1,2,3,4,5 };
    std::vector<double> y;
    double tolerance = 1e-12;
    for (double xi : x) {
        y.push_back(1 / (1 + xi * xi));
    }

    TridiagonalMatrix matrix(x.size());
    matrix.setValues_(x, y);

    std::vector<double> z = matrix.solveLU_(matrix.rightVector);
    double stepWidth = 1e-2;

    splineValues splineValues = calculateSplines_(x, y, z, stepWidth);

    std::vector<size_t> iteration;
    std::vector<double> omega;
    std::pair<size_t, std::vector<double>> itSol;
    for (double i = 0; i < 2; i += 1e-2) {
        omega.push_back(i);
        itSol = matrix.solveSOR_(matrix.rightVector, i, tolerance);
        iteration.push_back(itSol.first);
    }
    auto minIteraton = std::min_element(iteration.begin()+1, iteration.end());
    size_t indexMin = std::distance(iteration.begin(), minIteraton);
    std::cout << indexMin << std::endl;
    double omegaofMin = omega[indexMin];

    std::cout << "Minimum of Iteration: " << indexMin << " --> corresponding omega = " << omegaofMin << std::endl;
    plotResult2c(omega, iteration);
    plotResult2b(splineValues);
    */

    //-------- exercise 3 ----------------
    /*
    std::vector<std::vector<double>> data = readTxt2Matrix_("xyzm_dna.txt");

    std::vector<double> x = data[0];
    std::vector<double> y = data[1];
    std::vector<double> z = data[2];
    std::vector<double> m = data[3];
    size_t length = x.size();

    const int rCut = 5;
    const int k = 1;

    std::vector<std::vector<double>> hessianMatrix(length, std::vector<double>(length, 0.0));;

    for (size_t i = 0; i < length; ++i) {
        size_t neighbors = 0;

        for (size_t j = 0; j < length; ++j) {
            if (i != j) {
                double dx = x[i] - x[j];
                double dy = y[i] - y[j];
                double dz = z[i] - z[j];
                double rSquared = dx * dx + dy * dy + dz * dz;

                if (rSquared < (rCut * rCut)) {
                    ++neighbors;
                    hessianMatrix[i][j] = -k;
                }
            }

        }
        hessianMatrix[i][i] = neighbors * k;
    }

    std::vector<std::vector<double>> invSqrtmassMatrix(length, std::vector<double>(length, 0.0));
    for (size_t i = 0; i < length; ++i) {
        invSqrtmassMatrix[i][i] = 1/std::sqrt(m[i]);
    }

    std::vector<std::vector<double>> stiffnessMatrix = matrixMultiplication_(matrixMultiplication_(invSqrtmassMatrix, hessianMatrix) , invSqrtmassMatrix);
    std::vector<std::vector<double>> stiffnessMatrixDeflation(length, std::vector<double>(length,0.0));


    size_t maxIterations = 1000;
    std::vector<double> lambda(length);
    std::vector<std::vector<double>> eigenVectorMatrix(length, std::vector<double>(length, 0.0));
    for (size_t i = 0; i < 10; ++i) {
        eigenVectorMatrix[i] = powerMethod_(stiffnessMatrix, maxIterations);
        lambda[i] = eigenValues_(stiffnessMatrix, eigenVectorMatrix[i]);
        stiffnessMatrixDeflation = matrixMatrixSubtraction_(stiffnessMatrix, scalarMatrixMultiplication_(lambda[i], vectorVector2MatrixMultiplication_(eigenVectorMatrix[i], eigenVectorMatrix[i])));
        stiffnessMatrix = stiffnessMatrixDeflation;
    }
    plotResult3c(z, eigenVectorMatrix, lambda);
    */

    //-------- exercise 4 ----------------



}