#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <vectormath.h>


class StepTimer {
public:
    void startTimer();
    void stopStoreTimer();
    std::vector<float> getTimes();
private:
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    std::chrono::duration<float> duration;
    std::vector<float> times;
};


template<typename T>
void print(const T scalar) {
    std::cout << scalar << std::endl;
}

template <typename T>
void printVector(const std::vector<T>& vector) {
    std::cout << '[';
    for (size_t i = 0; i < vector.size(); ++i)
        std::cout << vector[i] << ' ';
    std::cout << ']' << std::endl;
}

template <typename T>
void printMatrix(const std::vector<std::vector<T>>& matrix) {
    for (const auto row : matrix) {
        for (T value : row) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
}

struct charge {
    double theta; //Angular position in spherical coordinates
    double phi; // Angular position in spherical coordubates
    double gradPhi; // Gradient of energy with respect to phi
    double gradTheta; // Gradient of energy with respect to theta
};

//----- Functions for exercise 1 

/*
Computes the potential energy for a given position in an anharmonic oscillator.
@param x The position at which the potential is evaluated.
@param m The mass of the particle.
@param lambda The anharmonic coefficient.
@param omega The angular frequency of the harmonic oscillator.
@return The potential energy at position x.
*/
double potential_(double x, double m, double lambda, double omega);

/*
Computes the orthodrome (great-circle distance) between two charges on a unit sphere.
@param q1 First charge.
@param q2 Second charge.
@return The orthodrome distance.
*/
double orthodrome_(const charge& q1, const charge& q2);

/*
Computes the Euclidean distance between two charges in 3D space.
@param q1 First charge.
@param q2 Second charge.
@return The Euclidean distance.
*/
double rij_(const charge& q1, const charge& q2);

/*
Calculates the total interaction energy of a system of charges.
@param charges A vector of charges.
@return Total energy of the system.
*/
double totalEnergy_(std::vector<charge> charges);

/*
Computes numerical gradients of energy with respect to charge positions (theta and phi).
@param charges A vector of charges.
@param stepSize Step size for finite difference approximation.
*/
void numericalGradient_(std::vector<charge>& charges, double stepSize);

/*
Performs gradient descent optimization to minimize the system's energy.
@param charges A vector of charges.
@param stepSize Step size for gradient descent updates.
@param maxIterations Maximum number of iterations.
*/
void gradientDecent_(std::vector<charge>& charges, double stepSize, double maxIterations);

//------ Functions for exercise 3

/*
Calculates the total force in a dimensionless system.
@param rNew Position vector in the dimensionless system.
@param R1New Position of the first heavy body in the dimensionless system.
@param R2New Position of the second heavy body in the dimensionless system.
@param omegaNew Angular velocity vector in the dimensionless system.
@param vNew Velocity vector in the dimensionless system.
@param mu Mass ratio parameter.
@return Total force vector in the dimensionless system.
*/
std::vector<double> totalForceDimLess_(std::vector<double> rNew, std::vector<double> R1New, std::vector<double> R2New, std::vector<double> omegaNew, std::vector<double> vNew, double mu);


std::vector<std::vector<double>> forwardEuler_(std::vector<double> rNew, std::vector<double> R1New, std::vector<double> R2New, std::vector<double> omegaNew, std::vector<double> vNew, double mu, double dt, int n);


std::vector<std::vector<double>> rungeKutta4_(std::vector<double> rNew, std::vector<double> R1New, std::vector<double> R2New, std::vector<double> omegaNew, std::vector<double> vNew, double mu, double dt, int n);

/*
Finds the Lagrange points for a given mass ratio (mu) in the circular restricted three-body problem.
@param mu The mass ratio of the two primary bodies.
@return A vector of vectors, where each inner vector represents the coordinates of a Lagrange point.
*/
std::vector<std::vector<double>> lagrangePointFinder_(double mu);

// ------ Functions for exercise 4

/*
Performs a single step of the Runge-Kutta method for solving ODEs.
@param x The current x-value.
@param y The current state vector.
@param h The step size.
@param func The derivative function defining the ODE system.
@return The updated state vector after the step.
*/
std::vector<double> rungeKutta_(double x, std::vector<double> y, double h, std::function<std::vector<double>(double, std::vector<double>)>func);

/*
Defines the Schrödinger equation as a system of ODEs.
@param x The spatial coordinate.
@param m The mass of the particle.
@param lambda The anharmonic potential parameter.
@param omega The angular frequency of the harmonic oscillator.
@param hbar The reduced Planck's constant.
@param E The energy level.
@param y The state vector [psi, psi'].
@return The derivative vector [psi', psi''].
*/
std::vector<double> schroedinger_(double x, double m, double lambda, double omega, double hbar, double E, std::vector<double>& y);

/*
Integrates the Schrödinger equation from x=0 to x=L and returns the wave function value at L.
@param E The energy level.
@param m The mass of the particle.
@param lambda The anharmonic potential parameter.
@param omega The angular frequency of the harmonic oscillator.
@param hbar The reduced Planck's constant.
@param L The integration limit.
@param h The step size for integration.
@param even If true, assumes an even initial wave function; otherwise, odd.
@return The value of the wave function at x=L.
*/
double psiL_(double E, double m, double lambda, double omega, double hbar, double L, double h, bool even);

/*
Finds energy brackets where the Schrödinger wave function changes sign.
@param m The mass of the particle.
@param lambda The anharmonic potential parameter.
@param omega The angular frequency of the harmonic oscillator.
@param hbar The reduced Planck's constant.
@param L The integration limit.
@param h The step size for integration.
@param N The number of desired brackets.
@param even If true, targets even wave functions; otherwise, odd.
@return A vector of energy intervals where a sign change occurs.
*/
std::vector<std::vector<double>> findBracket_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even);

/*
Finds eigenenergies of the Schrödinger equation using the bisection method.
@param m The mass of the particle.
@param lambda The anharmonic potential parameter.
@param omega The angular frequency of the harmonic oscillator.
@param hbar The reduced Planck's constant.
@param L The integration limit.
@param h The step size for integration.
@param N The number of eigenenergies to find.
@param even If true, targets even wave functions; otherwise, odd.
@param tolerance The tolerance for eigenenergy convergence.
@return A vector of eigenenergies.
*/
std::vector<double> eigenEnergyFinder_(double m, double lambda, double omega, double hbar, double L, double h, double N, bool even, double tolerance);
