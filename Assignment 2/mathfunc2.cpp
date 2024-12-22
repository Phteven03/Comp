#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>
#include <vectormath.h>

#include "mathfunc2.h"
#include "matplot/matplot.h"
#include "fftw3.h"
#include "matrixMath.h"


Timer::Timer(){
    start = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(0.0f);
    end = start;
}

Timer::~Timer(){
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    float ms = duration.count() * 1000.0f;
    std::cout << "Timer took " << ms << "ms " << std::endl;
}

void StepTimer::startTimer() {
    start = std::chrono::high_resolution_clock::now();
}

void StepTimer::stopStoreTimer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    float ms = duration.count() * 1000.0f;
    
    times.push_back(ms);
}

std::vector<float> StepTimer::getTimes() {
    return times;
}

std::vector<std::complex<double>> FFT_(const std::vector<double>& values, StepTimer* stepTimer) {

    if (stepTimer) {
        stepTimer->startTimer();
    }

    size_t N = values.size();
    std::vector<std::complex<double>> FFTResult(N);

    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2.0 + 1));
    double* in = (double*)fftw_malloc(sizeof(double) * N);
    std::copy(values.begin(), values.end(), in);

    fftw_plan plan = fftw_plan_dft_r2c_1d((int)N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    for (size_t i = 0; i < N / 2.0 + 1; ++i) {
        FFTResult[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    if (stepTimer) {
        stepTimer->stopStoreTimer();
    }

    return FFTResult;
}

std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double sampleRate, StepTimer* stepTimer) {
    std::vector<std::pair<double, double>> powerSpectrumVec;
    size_t N = values.size();

    if (stepTimer) {
        stepTimer->startTimer();
    }

    std::vector<std::complex<double>> FFTResult = FFT_(values, stepTimer);

    double freqResolution = sampleRate / N;

    for (size_t i = 0; i < N / 2 + 1; ++i) {

        double power = std::norm(FFTResult[i]);
        if (i > 0 && i < N / 2) {
            power *= 2.0;
        }
        double powerSpectrumValue = power * power / (N * N);

        double currentFreq = i * freqResolution;

        powerSpectrumVec.emplace_back(currentFreq, powerSpectrumValue);
    }

    if (stepTimer) {
        stepTimer->stopStoreTimer();
    }

    return powerSpectrumVec;
}

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth) {
    splineValues Sxi;

    size_t n = x.size();

    for (size_t i = 0; i < n-1; ++i) {
        double h = x[i + 1] - x[i];

        for (double xi = x[i]; xi <= x[i + 1]; xi += abs(x[i] - x[i + 1]) * stepWidth) {
            double splineValue = (z[i] / (6 * h)) * pow(x[i+1] - xi, 3) + (z[i+1] / (6 * h)) * pow(xi - x[i], 3) + ((y[i+1] / h) - (z[i+1] * h / 6)) * (xi - x[i]) + ((y[i] / h) - (z[i] * h / 6)) * (x[i+1] - xi);
            Sxi.splineValues.push_back(splineValue);
            Sxi.xValues.push_back(xi);
        }
    }
    return Sxi;
}

std::vector<double> powerMethod_(std::vector<std::vector<double>>& matrix, size_t& maxIterations) {
    size_t n = matrix[0].size();
    std::vector<double> guessVector(n, 1.0);

    for (size_t iterations = 0; iterations < maxIterations; ++iterations) {
        std::vector<double>nextVector = matrixVectorMuliplicaton_(matrix, guessVector);
        double norm = 0.0;
        double sum = 0.0;
        for (double val : nextVector) {
            sum += val * val;
        }
        norm = std::sqrt(sum);

        for (size_t i = 0; i < n; ++i) {
            nextVector[i] /= norm;
        }
        guessVector = nextVector;
    }
    return guessVector;
}

double eigenValues_(std::vector<std::vector<double>>& matrix, std::vector<double>& eigenVector) {
    double lambda = scalarProduct_(eigenVector, matrixVectorMuliplicaton_(matrix, eigenVector)) / scalarProduct_(eigenVector,eigenVector);
    return lambda;
}