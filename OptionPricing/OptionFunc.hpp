//	Header file for Option global functions

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "NumericalMethod.hpp"

using namespace std;

#ifndef OptionFunc_HPP
#define OptionFunc_HPP

double CallPayoff(double S, double K);	//	Call option payoff at maturity.
double PutPayoff(double S, double K);	//	Put option payoff at maturity.

double BSCallPrice(double S, double K, double T, double r, double sig, double q);	//	European call option pricing with Black-Scholes-Merton.
double BSPutPrice(double S, double K, double T, double r, double sig, double q);	//	European put option pricing with Black-Scholes-Merton.

double BSCallDelta(double S, double K, double T, double r, double sig, double q);	//	European call option Delta.
double BSPutDelta(double S, double K, double T, double r, double sig, double q);	//	European put option Delta.

double BSCallGamma(double S, double K, double T, double r, double sig, double q);	//	European call option Gamma.
double BSPutGamma(double S, double K, double T, double r, double sig, double q);	//	European put option Gamma.

double BSCallTheta(double S, double K, double T, double r, double sig, double q);	//	European call option Theta.
double BSPutTheta(double S, double K, double T, double r, double sig, double q);	//	European put option Theta.
											
double BSCallVega(double S, double K, double T, double r, double sig, double q);	//	European call option Vega.
double BSPutVega(double S, double K, double T, double r, double sig, double q);	//	European put option Vega.

vector<double> EuropeanBSPricer(string type, double S, double K, double T, double r, double sig, double q);	//	European Option Black-Scholes Pricer.

double BSCallImpVol(double V, double S, double K, double T, double r, double sig0, double q, double acc = 0.000001);	//	European call option Implied volatility with Black-Scholes-Merton.
double BSPutImpVol(double V, double S, double K, double T, double r, double sig0, double q, double acc = 0.000001);	//	European put option Implied volatility with Black-Scholes-Merton.

double BSDnOCallPrice(double S, double K, double T, double r, double sig, double q, double B);	//	Down-and-Out barrier call option pricing with closed formula.

struct Dividend
{
	double dividend;
	double time;
	string type;

	Dividend(double _dividend, double _time, string _type) : dividend(_dividend), time(_time), type(_type) {}
};

double DivToNon(vector<Dividend> q, double SDiv, double t, double T, double r);

double NonToDiv(vector<Dividend> q, double SNon, double t, double T, double r);

#endif