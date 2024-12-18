#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <iomanip>
#include <chrono>

#include "mathfunc2.h"
#include "matplot/matplot.h"
#include "fftw3.h"

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

std::vector<std::pair<double, double>> powerSpectrum_(const std::vector<double>& values, double frequency, StepTimer* stepTimer) {
    std::vector<std::pair<double, double>> powerSpectrumVec;
    size_t N = values.size();

    if (stepTimer) {
        stepTimer->startTimer();
    }

    std::vector<std::complex<double>> FFTResult = FFT_(values, stepTimer);

    double freqResolution = frequency / N;

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

TridiagonalMatrix createTridiagonalMatrix_(std::vector<double>& x, std::vector<double>& y) {
    size_t n = x.size();

    TridiagonalMatrix matrix;
    matrix.mid.resize(n-2);
    matrix.upper.resize(n - 3);
    matrix.rightvector.resize(n-2);

    std::vector<double> b(n - 1);
    std::vector<double> h(n - 1);

    for (size_t i = 0; i < n - 1 ; ++i) {
        h[i] = x[i + 1] - x[i];
        b[i] = (y[i + 1] - y[i]) / h[i];

        if (i < n - 3) {
            matrix.upper[i] = h[i];
        }
    }

    for (size_t i = 0; i < n-2; ++i) {
        matrix.mid[i] = 2 * (h[i] + h[i+1]);
        matrix.rightvector[i] = 6 * (b[i+1] - b[i]);
    }

    return matrix;
}

LUDecomposition LUD_(TridiagonalMatrix& matrix) {
    size_t n = matrix.mid.size();
    
    LUDecomposition LU;
    LU.UMid.resize(n);
    LU.UUpper.resize(n - 1);
    LU.LLower.resize(n - 1);

    LU.UMid[0] = matrix.mid[0];

    for (size_t i = 0; i < n - 1; ++i) {
        LU.UUpper[i] = matrix.upper[i];
    }

    for (size_t i = 1; i < n; ++i) {
        LU.LLower[i - 1] = matrix.upper[i - 1] / LU.UMid[i - 1];
        LU.UMid[i] = matrix.mid[i] - matrix.upper[i - 1] * LU.LLower[i - 1];
    }
    return LU;
}

std::vector<double> solveLU_(LUDecomposition& LU, std::vector<double>& u) {
    size_t n = u.size();
    std::vector<double> y(n);
    std::vector<double> z(n);

    y[0] = u[0];
    for(size_t i = 1; i < n; ++i) {
        y[i] = u[i] - LU.LLower[i-1] * y[i - 1];
    }

    z[n - 1] = y[n-1] / LU.UMid[n-1];
    for (int i = n - 2; i >= 0; --i) {
        z[i] = (y[i] - LU.UUpper[i] * z[i + 1]) / LU.UMid[i];
    }

    return z;
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

std::vector<double> matrixInversion_()

SOR calculateWithSOR_(TridiagonalMatrix& matrix, double& omega, double tolerance) {
    size_t n = matrix.rightvector.size();
    std::vector<double> guessSol(n, 0.0);
    std::vector<double> previousSol(n, 0.0);


}
