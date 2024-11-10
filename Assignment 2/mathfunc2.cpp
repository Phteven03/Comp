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
