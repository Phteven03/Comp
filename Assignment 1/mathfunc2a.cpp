#include <iostream>
#include <vector>
#include <random>
#include "mathfunc2a.h"

#include <iostream>
#include <vector>
#include <random>
#include "mathfunc2a.h"

double evaluatePoly_(std::vector<double> poly, double x) {
    double polyValue = 0.0;

    double power = 1.0;
    // polynomial evaluation using Horner's method
    for (double coeff : poly) {
        polyValue += coeff * power;
        power *= x;
    }

    return polyValue;
}

std::vector<double> polyDiff_(std::vector<double> polyFunction) {
    std::vector<double> polyDiff(polyFunction.size() - 1);

    // Compute the derivative of the polynomial
    for (int i = 1; i < polyFunction.size(); ++i) {
        polyDiff[i - 1] = i * polyFunction[i]; // Derivative formula
    }

    return polyDiff;
}

std::vector<double> bisectonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = bisection_(poly, significantIntervals);
    return convergenceSteps[0];
}

std::vector<double> newtonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = newtonRaphson_(poly, significantIntervals, 50);

    return convergenceSteps[0];
}

std::vector<double> polyRootBisection_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<double> roots;

    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = bisection_(poly, significantIntervals);

    // Extract the root for each convergence step
    for (const auto& steps : convergenceSteps) {
        if (!steps.empty()) {
            roots.push_back(steps.back());
        }
    }
    return roots;
}

std::vector<double> polyRootNewtonRaphson_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<double> roots;

    double n = 50;
    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = newtonRaphson_(poly, significantIntervals, n);

    // Extract the root for each convergence step
    for (const auto& steps : convergenceSteps) {
        if (!steps.empty()) {
            roots.push_back(steps.back());
        }
    }
    return roots;
}

std::vector<double> generateRandomNumbers(int numSamples) {
    std::vector<double> randomNumbers;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < numSamples; ++i) {
        randomNumbers.push_back(dis(gen));
    }

    return randomNumbers;
}

std::vector<std::vector<double>> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit) {
    std::vector<std::vector<double>> significantIntervalls;

    const double stepwidth = 1e-4;
    double previousValue = evaluatePoly_(poly, leftlimit);
    double previousX = leftlimit;

    // Look for intervals where the polynomial changes sign (potential roots)
    for (double currentX = previousX + stepwidth; currentX < rightlimit; currentX += stepwidth) {
        double currentValue = evaluatePoly_(poly, currentX);

        if ((previousValue * currentValue) < 0.0) {
            significantIntervalls.push_back({ previousX, currentX });
        }

        previousValue = currentValue;
        previousX = currentX;
    }

    return significantIntervalls;
}

std::vector<std::vector<double>> bisection_(std::vector<double> poly, std::vector<std::vector<double>> guessIntervalls) {
    std::vector<std::vector<double>> convergenceSteps;

    double tolerance = 1e-6;  // Tolerance for root approximation

    for (const auto& interval : guessIntervalls) {
        double leftLimit = interval[0];
        double rightLimit = interval[1];
        double midpoint;
        double leftValue = evaluatePoly_(poly, leftLimit);
        double rightValue = evaluatePoly_(poly, rightLimit);

        std::vector<double> intervalConvergence;

        // Bisection method to narrow down root intervals
        while ((rightLimit - leftLimit) / 2.0 > tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            double midpointValue = evaluatePoly_(poly, midpoint);
            intervalConvergence.push_back(midpoint);

            if (std::abs(midpointValue) < tolerance) {  // Root found within tolerance
                intervalConvergence.push_back(midpoint);
                break;
            }

            // Determine which side of the interval contains the root
            if (leftValue * midpointValue < 0) {
                rightLimit = midpoint;
                rightValue = midpointValue;
            }
            else {
                leftLimit = midpoint;
                leftValue = midpointValue;
            }
        }

        // Final midpoint when the interval is sufficiently small
        if ((rightLimit - leftLimit) / 2.0 <= tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            intervalConvergence.push_back(midpoint);
        }

        convergenceSteps.push_back(intervalConvergence);
    }

    return convergenceSteps;
}

std::vector<std::vector<double>> newtonRaphson_(std::vector<double> polyFunction, std::vector<std::vector<double>> guessIntervalls, int n) {
    std::vector<std::vector<double>> convergenceSteps;

    std::vector<double> f = polyFunction;
    std::vector<double> fPrime = polyDiff_(f);
    double tolerance = 1e-6;
    int maxIterations = n;

    for (const auto& interval : guessIntervalls) {
        double leftLimit = interval[0];
        double rightLimit = interval[1];
        double guess = (leftLimit + rightLimit) / 2.0;  // Initial guess is midpoint of the interval

        std::vector<double> intervalConvergence;

        // Newton-Raphson method for root-finding
        for (int i = 0; i < maxIterations; ++i) {
            double guessValue = evaluatePoly_(f, guess);
            double guessDiffValue = evaluatePoly_(fPrime, guess);
            double xNew = guess - guessValue / guessDiffValue;  // Newton-Raphson formula
            intervalConvergence.push_back(xNew);

            // If the difference is smaller than the tolerance, it converges
            if (std::abs(xNew - guess) < tolerance) {
                intervalConvergence.push_back(xNew);
                break;
            }
            guess = xNew;
        }

        convergenceSteps.push_back(intervalConvergence);
    }

    return convergenceSteps;
}


