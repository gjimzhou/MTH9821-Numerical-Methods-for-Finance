//	Header file for Trinomial Pricers

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "OptionFunc.hpp"

using namespace std;


#ifndef TrinomialTree_HPP
#define TrinomialTree_HPP

//	Trinomial Tree Pricer for European Options.
vector<double> EuropeanTrinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Trinomial Tree Pricer for American Options.
vector<double> AmericanTrinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Average Trinomial Tree Pricer for European Options.
vector<double> EuropeanATTPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Average Trinomial Tree Pricer for American Options.
vector<double> AmericanATTPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Trinomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanTBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Trinomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanTBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanTBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanTBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

//	Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Average Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanATTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Trinomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTBSVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);
//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTBSRVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N);

#endif


