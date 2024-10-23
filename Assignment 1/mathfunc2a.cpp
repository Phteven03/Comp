#include <iostream>
#include <vector>
#include <random>
#include "mathfunc2a.h"

double evaluatePoly_(std::vector<double> poly, double x) {
    double polyValue = 0.0;
    double power = 1.0;

    for (double coeff : poly) {
        polyValue += coeff * power;
        power *= x;
    }

    return polyValue;
}

std::vector<double> polyDiff_(std::vector<double> polyFunction) {
    std::vector<double> polyDiff(polyFunction.size() - 1);
    for (int i = 1; i < polyFunction.size(); ++i) {
        polyDiff[i - 1.0] = i * polyFunction[i];
    }
    return polyDiff;
}

std::vector<std::vector<double>> bracketing_(std::vector<double> poly, double leftlimit, double rightlimit) {
    const double stepwidth = 1e-4;
    std::vector<std::vector<double>> significantIntervalls;
    double previousValue = evaluatePoly_(poly, leftlimit);
    double previousX = leftlimit;
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
    double tolerance = 1e-6;

    for (const auto& interval : guessIntervalls) {
        double leftLimit = interval[0];
        double rightLimit = interval[1];
        double midpoint;
        double leftValue = evaluatePoly_(poly, leftLimit);
        double rightValue = evaluatePoly_(poly, rightLimit);

        std::vector<double> intervalConvergence;

        while ((rightLimit - leftLimit) / 2.0 > tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            double midpointValue = evaluatePoly_(poly, midpoint);
            intervalConvergence.push_back(midpoint);

            if (std::abs(midpointValue) < tolerance) {
                intervalConvergence.push_back(midpoint);
                break;
            }

            if (leftValue * midpointValue < 0) {
                rightLimit = midpoint;
                rightValue = midpointValue;
            }
            else {
                leftLimit = midpoint;
                leftValue = midpointValue;
            }
        }
        if ((rightLimit - leftLimit) / 2.0 <= tolerance) {
            midpoint = (leftLimit + rightLimit) / 2.0;
            intervalConvergence.push_back(midpoint);
        }
        convergenceSteps.push_back(intervalConvergence);
    }
    return convergenceSteps;
}

std::vector<double> bisectonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = bisection_(poly, { significantIntervals });
        return convergenceSteps[0];
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
        double guess = (leftLimit + rightLimit) / 2.0;

        std::vector<double> intervalConvergence; 
    
        for (int i = 0; i < maxIterations; ++i) {

            double guessValue = evaluatePoly_(f, guess);
            double guessDiffValue = evaluatePoly_(fPrime, guess);
            double xNew = guess - guessValue / guessDiffValue;
            intervalConvergence.push_back(xNew);

            if (abs(xNew - guess) < tolerance) {
                intervalConvergence.push_back(xNew);
                break;
            }
            guess = xNew;
        }

        convergenceSteps.push_back(intervalConvergence);
    }
    return convergenceSteps;
}

std::vector<double> newtonConvergence_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = newtonRaphson_(poly, { significantIntervals }, 100);

    return convergenceSteps[0];
}

std::vector<double> polyRootBisection_(std::vector<double> poly, double leftLimit, double rightLimit) {
    std::vector<double> roots;

    std::vector<std::vector<double>> significantIntervals = bracketing_(poly, leftLimit, rightLimit);
    std::vector<std::vector<double>> convergenceSteps = bisection_(poly, significantIntervals);

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

    for (const auto& steps : convergenceSteps) {
        if (!steps.empty()) {
            roots.push_back(steps.back());
        }
    }
    return roots;
}

/*std::vector<double> polyRootFinder_(std::vector<double> poly, double n, double leftlimit, double rightlimit) {
    std::vector<double> guess = bracketing_(poly, leftlimit, rightlimit);
    size_t numberOfRoots = guess.size();
    std::vector<double> polyRoot(numberOfRoots);
    polyRoot = newtonRaphson_(poly, guess, n);
    return polyRoot;
}*/

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