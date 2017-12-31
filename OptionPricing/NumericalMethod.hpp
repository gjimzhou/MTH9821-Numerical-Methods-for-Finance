//	Header file for Numerical methods functions

#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>

using namespace std;

#ifndef NumericalMethod_HPP
#define NumericalMethod_HPP

double pnorm(double x);	//	CDF for standard normal distribution.
double dnorm(double x);	//	PDF for standard normal distribution.
double qnorm(double x);	//	QF for standard normal distribution.

long int SeedInit(long int seed = 0);	//	Initialize random seed.

//	Linear Congruential Generator for standard uniform distribution random variable.
vector<double> LinearCongruentialGenerator(long int N, long int seed = 0);

// Inverse Transform Method for standard normal distribution random variable.
vector<double> InverseTransformMethod(long int N, long int seed = 0);

// Acceptance-Rejection Method for standard normal distribution random variable.
vector<double> AcceptanceRejectionMethod(long int N, long int seed = 0);

// Box-Muller Method for standard normal distribution random variable.
vector<double> BoxMullerMethod(long int N, long int seed = 0);

//	Generator for standard normal distribution random variables.
vector<double> rnorm(long int N, string Generator = "ITM", long int seed = 0);

#endif
