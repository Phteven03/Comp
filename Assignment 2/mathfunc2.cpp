#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>
#include <vectormath.h>
#include <algorithm>

#include "mathfunc2.h"
#include "matplot/matplot.h"
#include "fftw3.h"
#include "matrixMath.h"


Timer::Timer() {
    start = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(0.0f);
    end = start;
}

Timer::~Timer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;

    float ms = duration.count() * 1000.0f;
    std::cout << "Timer took " << ms << "ms " << std::endl;
}

void StepTimer::startTimer() {
    start = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(0.0f);
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

std::vector<std::complex<double>> FFT_(const std::vector<double>& values, int signExp, StepTimer* stepTimer) {
    if (stepTimer) {
        stepTimer->startTimer();
    }

    size_t N = values.size();
    std::vector<std::complex<double>> FFTResult(N);

    fftw_complex* in = fftw_alloc_complex(N);
    fftw_complex* out = fftw_alloc_complex(N);

    for (int i = 0; i < N; ++i)
    {
        in[i][0] = values[i];
        in[i][1] = 0.0;
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, signExp, FFTW_ESTIMATE);
    fftw_execute(plan);

    for (size_t i = 0; i < N; ++i) {
        FFTResult[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();

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

    std::vector<std::complex<double>> FFTResult = FFT_(values, -1, stepTimer);

    double freqResolution = sampleRate / N;

    for (size_t i = 0; i < N / 2 + 1; ++i) {

        double power = std::norm(FFTResult[i]);
        //if (i > 0 && i < N / 2) {
        //    power *= 2.0;
        //}
        double powerSpectrumValue = power * power / (N * N);

        double currentFreq = i * freqResolution;

        powerSpectrumVec.emplace_back(currentFreq, power);
    }

    if (stepTimer) {
        stepTimer->stopStoreTimer();
    }

    return powerSpectrumVec;
}

std::vector<double> bubbleSort_(std::vector<double> vector) {
    size_t n = vector.size();
    for (size_t j = n; j > 1; --j) {
        for (size_t i = 0; i < n - 1; ++i) {
            if (vector[i] < vector[i + 1]) {
                std::swap(vector[i], vector[i + 1]);
            }
        }
    }
    return vector;
}

std::vector<double> maxFinder_(std::vector<double> vector) {
    size_t n = vector.size();
    for (size_t i = 0; i < n - 1; ++i) {
        if (vector[i] < vector[i + 1]) {
            vector[i] = 0.0;
        }
        if (vector[i] > vector[i + 1]) {
            vector[i + 1] = 0.0;
        }
    }
    for (size_t i = n; i > 0; --i) {
        vector[i] = vector[i] < 1600 ? vector[i] = 0.0 : vector[i];
    }
    vector.erase(std::remove(vector.begin(), vector.end(), 0.0), vector.end());

    return vector;
}

splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth) {
    splineValues Sxi;

    size_t n = x.size();

    for (size_t i = 0; i < n - 1; ++i) {
        double h = x[i + 1] - x[i];

        for (double xi = x[i]; xi <= x[i + 1]; xi += abs(x[i] - x[i + 1]) * stepWidth) {
            double splineValue = (z[i] / (6 * h)) * pow(x[i + 1] - xi, 3) + (z[i + 1] / (6 * h)) * pow(xi - x[i], 3) + ((y[i + 1] / h) - (z[i + 1] * h / 6)) * (xi - x[i]) + ((y[i] / h) - (z[i] * h / 6)) * (x[i + 1] - xi);
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
    double lambda = scalarProduct_(eigenVector, matrixVectorMuliplicaton_(matrix, eigenVector)) / scalarProduct_(eigenVector, eigenVector);
    return lambda;
}