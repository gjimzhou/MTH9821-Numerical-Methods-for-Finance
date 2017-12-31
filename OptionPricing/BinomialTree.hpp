//	Header file for Binomial Pricers

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "OptionFunc.hpp"

using namespace std;


#ifndef BinomialTree_HPP
#define BinomialTree_HPP

//	Binomial Tree Pricer for European Options.
vector<double> EuropeanBinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Binomial Tree Pricer for American Options.
vector<double> AmericanBinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Average Binomial Tree Pricer for European Options.
vector<double> EuropeanABTPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Average Binomial Tree Pricer for American Options.
vector<double> AmericanABTPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Binomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanBBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Binomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanBBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanBBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanBBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Average Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanABTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Binomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBBSVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBBSRVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Implied Volatility for American Options using Binomial Tree Pricer.
double AmericanBTImpVol(string type, double Vm, double S, double K, double T, double r, double q, int N, double sig0 = 0.1, double sig1 = 0.5, double thres = 0.00001);


//	Binomial Tree Pricer for European Options with Discrete Dividends.
vector<double> EuropeanDividendBinomialTreePricer(string type, double S, double K, double T, double r, double sig, vector<Dividend> q, int N);
//	Binomial Tree Pricer for American Options with Discrete Dividends.
vector<double> AmericanDividendBinomialTreePricer(string type, double S, double K, double T, double r, double sig, vector<Dividend> q, int N);


#endif

