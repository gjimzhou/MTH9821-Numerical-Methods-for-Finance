//	Header file for Option class

#include <iostream>
#include <string>
#include <cmath>

#include "OptionFunc.hpp"
#include "BinomialTree.hpp"
#include "TrinomialTree.hpp"
#include "MonteCarlo.hpp"

using namespace std;

#ifndef Option_HPP
#define Option_HPP

// Parent Option
class Option
{
protected:
	//	Option data (K, T, r, sig, q).
	double K;
	double T;
	double r;
	double sig;
	double q;
	string type = "C";	//	Option type (C or P).
	string name = "APPL";	//	Asset name.

public:
	Option();	//	Default constructor.
	Option(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type = "C", string new_name = "APPL");	//	Constructor.
	Option(const Option& o);	//	Copy constructor.
	virtual ~Option();	//	Destructor.

	Option& operator = (const Option& o);	//	Assignment Operator.

	double Strike();	//	Get K.
	void Strike(double new_K);	//	Set K.
	double Maturity();	//	Get T.
	void Maturity(double new_T);	//	Set T.
	double InterestRate();	//	Get r.
	void InterestRate(double new_r);	//	Set r.
	virtual double Volatility();	//	Get sig.
	virtual void Volatility(double new_sig);	//	Set sig.
	virtual double Dividend();	//	Get q.
	virtual void Dividend(double new_q);	//	Set q.

	string Name();	//	Get name.
	void Name(string new_name);	//	Set name.

	string Type();	//	Get type.
	void Type(string new_type);	//	Set type.
	void Toggle();	//	Switch type.

};


//	European Option
class EuropeanOption : public Option
{
public:
	EuropeanOption();	//	Default constructor.
	EuropeanOption(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type = "C", string new_name = "APPL");	//	Constructor.
	EuropeanOption(const EuropeanOption& o);	//	Copy constructor.
	virtual ~EuropeanOption();	//	Destructor.

	EuropeanOption& operator = (const EuropeanOption& o);	//	Assignment Operator.

	double BSPrice(double S);	//	Calculate price with Black-Scholes-Merton.
	double BSDelta(double S);	//	Calculate Delta with Black-Scholes-Merton.
	double BSGamma(double S);	//	Calculate Gamma with Black-Scholes-Merton.
	double BSTheta(double S);	//	Calculate Theta with Black-Scholes-Merton.
	vector<double> BSPricer(double S);	//	European Option Black-Scholes Pricer.

	double BSImpVol(double V, double S, double acc = 0.000001);	//	Calculate Implied Volatility with Black-Scholes-Merton.

	vector<double> BinomialTreePricer(double S, int N);	//	Binomial Tree Pricer for European Options.
	vector<double> ABTPricer(double S, int N);	//	Average Binomial Tree Pricer for European Options.
	vector<double> BBSPricer(double S, int N);	//	Binomial Tree Black-Scholes Pricer for European Options.
	vector<double> BBSRPricer(double S, int N);	//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.

	vector<double> TrinomialTreePricer(double S, int N);	//	Trinomial Tree Pricer for European Options.
	vector<double> ATTPricer(double S, int N);	//	Average Trinomial Tree Pricer for European Options.
	vector<double> TBSPricer(double S, int N);	//	Trinomial Tree Black-Scholes Pricer for European Options.
	vector<double> TBSRPricer(double S, int N);	//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.

	vector<double> MonteCarloPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.
	
	double MCCVPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Control Variate Technique.
	double MCAVPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Antithetic Variablews.
	double MCMMPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Moment Matching.
	double MCCVMMPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Simultaneous Moment Matching and Control Variate Technique.

};


//	American Option
class AmericanOption : public Option
{
public:
	AmericanOption();	//	Default constructor.
	AmericanOption(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type = "C", string new_name = "APPL");	//	Constructor.
	AmericanOption(const AmericanOption& o);	//	Copy constructor.
	virtual ~AmericanOption();	//	Destructor.

	AmericanOption& operator = (const AmericanOption& o);	//	Assignment Operator.

	vector<double> BinomialTreePricer(double S, int N);	//	Binomial Tree Pricer for American Options.
	vector<double> ABTPricer(double S, int N);	//	Average Binomial Tree Pricer for American Options.
	vector<double> BBSPricer(double S, int N);	//	Binomial Tree Black-Scholes Pricer for American Options.
	vector<double> BBSRPricer(double S, int N);	//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
	
	vector<double> BTVRPricer(double S, int N);	//	Binomial Tree Pricer for American Options.	//	With variance reduction.
	vector<double> ABTVRPricer(double S, int N);	//	Average Binomial Tree Pricer for American Options.	//	With variance reduction.
	vector<double> BBSVRPricer(double S, int N);	//	Binomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
	vector<double> BBSRVRPricer(double S, int N);	//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.

	vector<double> TrinomialTreePricer(double S, int N);	//	Trinomial Tree Pricer for American Options.
	vector<double> ATTPricer(double S, int N);	//	Average Trinomial Tree Pricer for American Options.
	vector<double> TBSPricer(double S, int N);	//	Trinomial Tree Black-Scholes Pricer for American Options.
	vector<double> TBSRPricer(double S, int N);	//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.

	vector<double> TTVRPricer(double S, int N);	//	Trinomial Tree Pricer for American Options.	//	With variance reduction.
	vector<double> ATTVRPricer(double S, int N);	//	Average Trinomial Tree Pricer for American Options.	//	With variance reduction.
	vector<double> TBSVRPricer(double S, int N);	//	Trinomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
	vector<double> TBSRVRPricer(double S, int N);	//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.

};


//	Barrier Option
class BarrierOption : public Option
{
protected:
	//	Barrier (B, Btype).
	double B;
	string Btype = "DnO";

public:
	BarrierOption();	//	Default constructor.
	BarrierOption(double new_K, double new_T, double new_r, double new_sig, double new_q, double new_B,  string new_type = "C", string new_Btype = "DnO", string new_name = "APPL");	//	Constructor.
	BarrierOption(const BarrierOption& o);	//	Copy constructor.
	virtual ~BarrierOption();	//	Destructor.

	BarrierOption& operator = (const BarrierOption& o);	//	Assignment Operator.

	double Barrier();	//	Get B.
	void Barrier(double new_B);	//	Set B.

	string BarrierType();	//	Get Btype.
	void BarrierType(string new_Btype);	//	Set Btype.

	double BSPrice(double S);	//	Calculate price with Black-Scholes-Merton.

	double MonteCarloPricer(double S, long int n, long int m, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for Barrier Options.
	double MonteCarloPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for Barrier Options.	//	With optimal path number and time interval number.

};


//	Basket Option
class BasketOption : public Option
{
protected:
	//	Basket Option data (sig2, q2, rho).
	double sig2;
	double q2;
	double rho;
	string type = "C";	//	Option type (C or P).
	string name = "APPL";	//	Asset name.

public:
	BasketOption();	//	Default constructor.
	BasketOption(double new_K, double new_T, double new_r, double new_sig1, double new_sig2, double new_q1, double new_q2, double new_rho, string new_type = "C", string new_name = "APPL");	//	Constructor.
	BasketOption(const BasketOption& o);	//	Copy constructor.
	virtual ~BasketOption();	//	Destructor.

	BasketOption& operator = (const BasketOption& o);	//	Assignment Operator.

	double Volatility(bool n);	//	Get siq.
	void Volatility(double new_sig, bool n);	//	Set sig.

	double Dividend(bool n);	//	Get q.
	void Dividend(double new_q, bool n);	//	Set q.

	double Rho();	//	Get rho.
	void Rho(double new_rho);	//	Set rho.

	virtual double MonteCarloPricer(double S1, double S2, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for Basket Options.

};


//	Lookback Basket Option
class LookbackBasketOption : public BasketOption
{
public:
	LookbackBasketOption();	//	Default constructor.
	LookbackBasketOption(double new_K, double new_T, double new_r, double new_sig1, double new_sig2, double new_q1, double new_q2, double new_rho, string new_type = "C", string new_name = "APPL");	//	Constructor.
	LookbackBasketOption(const LookbackBasketOption& o);	//	Copy constructor.
	virtual ~LookbackBasketOption();	//	Destructor.

	LookbackBasketOption& operator = (const LookbackBasketOption& o);	//	Assignment Operator.

	double MonteCarloPricer(double S1, double S2, long int n, long int m, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for Lookback Basket Options.
	double MonteCarloPricer(double S1, double S2, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for Lookback Basket Options.	//	With optimal path number and time interval number.

};


//	European Option with stock price follows Heston Model
class EuropeanHestonOption : public EuropeanOption
{
protected:
	//	Heston model coefficient (lambda, sigbar, eta, rho).
	double lambda;
	double sigbar;
	double eta;
	double rho;

public:
	EuropeanHestonOption();	//	Default constructor.
	EuropeanHestonOption(double new_K, double new_T, double new_r, double new_sig, double new_lambda, double new_sigbar, double new_eta, double new_rho, double new_q, string new_type = "C", string new_name = "APPL");	//	Constructor.
	EuropeanHestonOption(const EuropeanHestonOption& o);	//	Copy constructor.
	virtual ~EuropeanHestonOption();	//	Destructor.

	EuropeanHestonOption& operator = (const EuropeanHestonOption& o);	//	Assignment Operator.

	double Lambda();	//	Get lambda.
	void Lambda(double new_lambda);	//	Set lambda.

	double MeanVolatility();	//	Get sigbar.
	void MeanVolatility(double new_sigbar);	//	Set sigbar.

	double Eta();	//	Get eta.
	void Eta(double new_eta);	//	Set eta.

	double Rho();	//	Get rho.
	void Rho(double new_rho);	//	Set rho.

	double MonteCarloPricer(double S, long int n, long int m, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.
	double MonteCarloPricer(double S, long int N, string Generator = "ITM", long int seed = 0);	//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.	//	With optimal path number and time interval number.

};


#endif