#pragma once
#include <iostream>
#include <vector>

double evaluatePoly_(std::vector<double> poly, double x);

std::vector<double> polyDiff_(std::vector<double> polyFunction);

std::vector<std::vector<double>> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit);

std::vector<double> bisection_(std::vector<double> poly, std::vector<std::vector<double>> guessIntervalls);

std::vector<double> newtonRaphson_(std::vector<double> polyFunction, std::vector<std::vector<double>> guessIntervalls, int n);

std::vector<double> polyRootBisection_(std::vector<double> poly, double leftLimit, double rightLimit);

std::vector<double> polyRootNewtonRaphson_(std::vector<double> poly, double leftLimit, double rightLimit);

//std::vector<double> polyRootFinder_(std::vector<double> poly, double n, double leftlimit, double rightlimit); 