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


//	Transpose a matrix.
vector<vector<double>> transpose(vector<vector<double>> A);

//	Forward Substitution.
vector<double> ForwardSubstitution(vector<vector<double>> L, vector<double> b);

//	Backward Substitution.
vector<double> BackwardSubstitution(vector<vector<double>> U, vector<double> b);

//	Forward Substitution for bidiagonal matrix.
vector<double> ForwardSubstitutionBidiag(vector<vector<double>> L, vector<double> b);

//	Backward Substitution for bidiagonal matrix.
vector<double> BackwardSubstitutionBidiag(vector<vector<double>> U, vector<double> b);


//	LU Decomposition without pivoting.
vector<vector<vector<double>>> LUNoPivoting(vector<vector<double>> A);

//	LU Decomposition without pivoting for tridiagonal matrix.
vector<vector<vector<double>>> LUNoPivotingTridiag(vector<vector<double>> A);

//	Linear Solver using LU Decomposition without pivoting.
vector<double> LinearSolverLUNoPivoting(vector<vector<double>> A, vector<double> b);

//	Linear Solver using LU Decomposition without pivoting for tridiagonal matrix.
vector<double> LinearSolverLUNoPivotingTridiag(vector<vector<double>> A, vector<double> b);


//	Cholesky Decomposition.
vector<vector<double>> Cholesky(vector<vector<double>> A);

//	Linear Solver using Cholesky Decomposition.
vector<double> LinearSolverCholesky(vector<vector<double>> A, vector<double> b);


#endif
