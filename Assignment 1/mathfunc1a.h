#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>

struct Timer {
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;

	Timer();

	~Timer();
};
template <typename T>
void printVector(std::vector<T>& vector) {
    std::cout << '[';
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << vector[i] << ' ';
    std::cout << ']' << std::endl;
}

double legendrePn_(double x, int n);

double legendrePnDiff_(double x, int n);

std::vector<double> newtonRaphson_(std::vector<double> guesses, double n);

std::vector<double> bracketing_(double leftlimit, double rightlimit, double n);

std::vector<double> gaussLegendreWeight_(std::vector<double> roots, double n);

std::vector<double> legendreRootFinder_(double leftlimit, double rightlimit, double n);