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
    matrix.mid.resize(n-1);
    matrix.upper.resize(n - 2);
    matrix.rightvector.resize(n-1);

    std::vector<double> b(n - 1);
    std::vector<double> h(n - 1);

    for (size_t i = 0; i < n - 1 ; ++i) {
        h[i] = x[i + 1] - x[i];
        b[i] = (y[i + 1] - y[i]) / h[i];

        if (i < n - 2) {
            matrix.upper[i] = h[i];
        }
    }

    for (size_t i = 0; i <= n-1; ++i) {
        matrix.mid[i + 1] = 2 * (h[i] + h[i+1]);
        matrix.rightvector[i] = 6 * (b[i+1] - b[i]);
    }

    return matrix;
}

//LUDecomposition LUD_(TridiagonalMatrix& matrix) {
//    size_t n = matrix.mid.size();
//    
//    LUDecomposition LU;
//    LU.UMid.resize(n);
//    LU.UUpper.resize(n - 1);
//    LU.LLower.resize(n - 1);
//
//    LU.UMid[0] = matrix.mid[0];
//
//    for (size_t i = 0; i < n - 1; ++i) {
//        LU.UUpper[i] = matrix.upper[i];
//    }
//
//    for (size_t i = 1; i < n; ++i) {
//        LU.LLower[i - 1] = matrix.upper[i - 1] / LU.UMid[i - 1];
//        LU.UMid[i] = matrix.mid[i] - matrix.upper[i - 1] * LU.LLower[i - 1];
//    }
//    return LU;
//}

//std::vector<double> solveLU_(LUDecomposition& LU, std::vector<double>& u) {
//    size_t n = u.size();
//    std::vector<double> y(n);
//    std::vector<double> z(n);
//
//    y[0] = u[0] / LU.UMid[0];
//    for (size_t i = 1; i < n; ++i) {
//        y[i] = (u[i] - LU.LLower[i - 1] * y[i - 1]) / LU.UMid[i];
//    }
//
//    z[n - 1] = y[n - 1];
//    for (int i = n - 2; i >= 0; --i) {
//        z[i] = y[i] - LU.UUpper[i] * z[i + 1] / LU.UMid[i];
//    }
//    return z;
//}
//
//double evaluateSplineSegment_(double x, double x_i, double x_ip1, double y_i, double y_ip1, double z_i, double z_ip1, double h) {
//    double term1 = (z_i / (6 * h)) * pow(x_ip1 - x, 3);
//    double term2 = (z_ip1 / (6 * h)) * pow(x - x_i, 3);
//    double term3 = ((y_ip1 / h) - (z_ip1 * h / 6)) * (x - x_i);
//    double term4 = ((y_i / h) - (z_i * h / 6)) * (x_ip1 - x);
//    return term1 + term2 + term3 + term4;
//}
//
//splineValues calculateSplines_(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double stepWidth) {
//    splineValues Sxi;
//
//    size_t n = x.size() - 1;
//
//    for (size_t i = 0; i < n; ++i) {
//        double h = x[i + 1] - x[i];
//
//        for (double xi = x[i]; xi <= x[i + 1]; xi += stepWidth) {
//            double splineValue = evaluateSplineSegment_(xi, x[i], x[i + 1], y[i], y[i + 1], z[i], z[i + 1], h);
//            Sxi.splineValues.push_back(splineValue);
//            Sxi.xValues.push_back(xi);
//        }
//    }
//    return Sxi;
//}
