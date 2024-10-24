#pragma once
#include <iostream>
#include <vector>
#include <random>

/*
  Evaluates a polynomial at a given point x
  @param poly Coefficients of the polynomial
  @param x Point at which to evaluate the polynomial
  @return The value of the polynomial at x
*/
double evaluatePoly_(std::vector<double> poly, double x);

/*
  Computes the derivative of a polynomial
  @param polyFunction Coefficients of the polynomial
  @return Coefficients of the derived polynomial
*/
std::vector<double> polyDiff_(std::vector<double> polyFunction);

/*
  Finds intervals where the polynomial changes sign (potential roots)
  @param poly Coefficients of the polynomial
  @param leftlimit Left boundary of the search range
  @param rightlimit Right boundary of the search range
  @return Intervals where roots may be located
*/
std::vector<std::vector<double>> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit);

/*
  Performs the bisection method to narrow down root intervals
  @param poly Coefficients of the polynomial
  @param guessIntervalls Initial intervals where roots are suspected
  @return Refined intervals after applying the bisection method
*/
std::vector<std::vector<double>> bisection_(std::vector<double> poly, std::vector<std::vector<double>> guessIntervalls);

/*
  Applies the bisection method to converge on a root
  @param poly Coefficients of the polynomial
  @param leftLimit Left boundary of the interval
  @param rightLimit Right boundary of the interval
  @return Approximate root of the polynomial
*/
std::vector<double> bisectonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit);

/*
  Performs the Newton-Raphson method to refine root intervals
  @param polyFunction Coefficients of the polynomial
  @param guessIntervalls Initial intervals where roots are suspected
  @param n Number of iterations for the method
  @return Refined intervals after applying Newton-Raphson
*/
std::vector<std::vector<double>> newtonRaphson_(std::vector<double> polyFunction, std::vector<std::vector<double>> guessIntervalls, int n);

/*
  Applies the Newton-Raphson method to converge on a root
  @param poly Coefficients of the polynomial
  @param leftLimit Left boundary of the interval
  @param rightLimit Right boundary of the interval
  @return Approximate root of the polynomial
*/
std::vector<double> newtonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit);

/*
  Combines bisection methods to find a root of the polynomial
  @param poly Coefficients of the polynomial
  @param leftLimit Left boundary of the interval
  @param rightLimit Right boundary of the interval
  @return Approximate root of the polynomial
*/
std::vector<double> polyRootBisection_(std::vector<double> poly, double leftLimit, double rightLimit);

/*
  Combines Newton-Raphson methods to find a root of the polynomial
  @param poly Coefficients of the polynomial
  @param leftLimit Left boundary of the interval
  @param rightLimit Right boundary of the interval
  @return Approximate root of the polynomial
*/
std::vector<double> polyRootNewtonRaphson_(std::vector<double> poly, double leftLimit, double rightLimit);

/*
  Generates a list of random numbers
  @param numSamples Number of equaly distant random numbers to generate
  @return A vector of random numbers
*/
std::vector<double> generateRandomNumbers(int numSamples);
