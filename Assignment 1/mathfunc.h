#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

double factorial_(double x);

double logFactorial_(int n);

double evaluatePoly_(std::vector<double> poly, double x);

std::vector<double> legendrePn_(int l);

std::vector<double> polyDiff_(std::vector<double> polyFunction);

std::vector<double> newtonRaphson_(std::vector<double> polyFunction, std::vector<double> guesses, int n);

std::vector<double> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit);

std::vector<double> gaussLegendreWeight_(std::vector<double> legendrePoly, std::vector<double> roots);

std::vector<double> polyRootFinder_(std::vector<double> poly, double n, double leftlimit, double rightlimit);