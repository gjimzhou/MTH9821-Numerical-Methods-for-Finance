//	Header file for Monte-Carlo Pricers

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "NumericalMethod.hpp"
#include "OptionFunc.hpp"

using namespace std;


#ifndef MonteCarlo_HPP
#define MonteCarlo_HPP

//	Monte-Carlo Pricer for European Options.
vector<double> EuropeanMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for Barrier Options.
double BarrierMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, double B, string Btype, long int n, long int m, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for Barrier Options.	//	With optimal path number and time interval number.
double BarrierMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, double B, string Btype, long int N, string Generator = "ITM", long int seed = 0);


//	Monte-Carlo Pricer for European Options.	//	Control Variate Technique.
double EuropeanMCCVPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for European Options.	//	Antithetic Variablews.
double EuropeanMCAVPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for European Options.	//	Moment Matching.
double EuropeanMCMMPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for European Options.	//	Simultaneous Moment Matching and Control Variate Technique.
double EuropeanMCCVMMPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator = "ITM", long int seed = 0);


//	Monte-Carlo Pricer for European Basket Options.
double BasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int N, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for Lookback Basket Options.
double LookbackBasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int n, long int m, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for Lookback Basket Options.	//	With optimal path number and time interval number.
double LookbackBasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int N, string Generator = "ITM", long int seed = 0);


//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.
double EuropeanHestonMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double lambda, double sigbar, double eta, double rho, double q, long int n, long int m, string Generator = "ITM", long int seed = 0);

//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.	//	With optimal path number and time interval number.
double EuropeanHestonMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double lambda, double sigbar, double eta, double rho, double q, long int N, string Generator = "ITM", long int seed = 0);


#endif
