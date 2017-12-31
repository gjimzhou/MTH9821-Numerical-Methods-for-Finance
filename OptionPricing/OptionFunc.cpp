//	Source file for Option global functions

#include "OptionFunc.hpp"

using namespace std;


//	Call option payoff at maturity.
double CallPayoff(double S, double K)
{
	double payoff = S - K;
	return payoff > 0 ? payoff : 0;
}

//	Put option payoff at maturity.
double PutPayoff(double S, double K)
{
	double payoff = K - S;
	return payoff > 0 ? payoff : 0;
}


//	European call option pricing with Black-Scholes-Merton.
double BSCallPrice(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);

	return (S * exp(-q * T) * pnorm(d1)) - (K * exp(-r * T) * pnorm(d2));
}

//	European put option pricing with Black-Scholes-Merton.
double BSPutPrice(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);

	return (K * exp(-r * T) * pnorm(-d2)) - (S * exp(-q * T) * pnorm(-d1));
}


//	European call option Delta.
double BSCallDelta(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return exp(-q * T) * pnorm(d1);
}

//	European put option Delta.
double BSPutDelta(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return -exp(-q * T) * pnorm(-d1);
}

//	European call option Gamma.
double BSCallGamma(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return (exp(-q * T) * dnorm(d1)) / (S * sig * sqrt(T));
}

//	European put option Gamma.
double BSPutGamma(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return (exp(-q * T) * dnorm(d1)) / (S * sig * sqrt(T));
}

//	European call option Theta.
double BSCallTheta(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);
	return -(S * sig * exp(-q * T) * dnorm(d1)) / (2 * sqrt(T)) + q * S * exp(-q * T) * pnorm(d1) - r * K * exp(-r * T) * pnorm(d2);
}

//	European put option Theta.
double BSPutTheta(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	double d2 = d1 - sig * sqrt(T);
	return -(S * sig * exp(-q * T) * dnorm(d1)) / (2 * sqrt(T)) - q * S * exp(-q * T) * pnorm(-d1) + r * K * exp(-r * T) * pnorm(-d2);
}

//	European call option Vega.
double BSCallVega(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return S * exp(-q * T) * sqrt(T) * dnorm(d1);
}

//	European put option Vega.
double BSPutVega(double S, double K, double T, double r, double sig, double q)
{
	double d1 = (log(S / K) + (r - q + (sig * sig) * 0.5) * T) / (sig * sqrt(T));
	return S * exp(-q * T) * sqrt(T) * dnorm(d1);
}

//	European Option Black-Scholes Pricer.
vector<double> EuropeanBSPricer(string type, double S, double K, double T, double r, double sig, double q)
{
	vector<double> result;
	double V, Delta, Gamma, Theta, Vega;

	if (type == "C")
	{
		V = BSCallPrice(S, K, T, r, sig, q);
		Delta = BSCallDelta(S, K, T, r, sig, q);
		Gamma = BSCallGamma(S, K, T, r, sig, q);
		Theta = BSCallTheta(S, K, T, r, sig, q);
		Vega = BSCallVega(S, K, T, r, sig, q);
	}
	else
	{
		V = BSPutPrice(S, K, T, r, sig, q);
		Delta = BSPutDelta(S, K, T, r, sig, q);
		Gamma = BSPutGamma(S, K, T, r, sig, q);
		Theta = BSPutTheta(S, K, T, r, sig, q);
		Vega = BSPutVega(S, K, T, r, sig, q);
	}

	result.push_back(V);
	result.push_back(Delta);
	result.push_back(Gamma);
	result.push_back(Theta);
	result.push_back(Vega);

	return result;
}


//	European call option Implied volatility with Black-Scholes-Merton.
double BSCallImpVol(double V, double S, double K, double T, double r, double sig0, double q, double acc)
{
	double sig = sig0 + 2 * acc;

	while (abs(sig - sig0) > acc)
	{
		sig = sig0;
		sig0 = sig0 - (BSCallPrice(S, K, T, r, sig0, q) - V) / BSCallVega(S, K, T, r, sig0, q);
	}

	return sig0;
}

//	European put option Implied volatility with Black-Scholes-Merton.
double BSPutImpVol(double V, double S, double K, double T, double r, double sig0, double q, double acc)
{
	double sig = sig0 + 2 * acc;

	while (abs(sig - sig0) > acc)
	{
		sig = sig0;
		sig0 = sig0 - (BSPutPrice(S, K, T, r, sig0, q) - V) / BSPutVega(S, K, T, r, sig0, q);
	}

	return sig0;
}


//	Down-and-Out barrier call option pricing with closed formula.
double BSDnOCallPrice(double S, double K, double T, double r, double sig, double q, double B)
{
	double a = (r - q) / (sig * sig) - 0.5; 
	return BSCallPrice(S, K, T, r, sig, q) - pow(B / S, 2 * a) * BSCallPrice(B * B / S, K, T, r, sig, q);
}


double DivToNon(vector<Dividend> q, double SDiv, double t, double T, double r)
{
	for (int i = 0; i < q.size(); i++)
	{
		if (t <= q[i].time && T >= q[i].time)
		{
			if (q[i].type == "p")
			{
				SDiv *= q[i].dividend;
			}
			else
			{
				SDiv -= q[i].dividend * exp(-r * (q[i].time - t));
			}
		}
	}

	return SDiv;
}

double NonToDiv(vector<Dividend> q, double SNon, double t, double T, double r)
{
	for (int i = q.size() - 1; i >= 0; i--)
	{
		if (t <= q[i].time && T >= q[i].time)
		{
			if (q[i].type == "p")
			{
				SNon /= q[i].dividend;
			}
			else
			{
				SNon += q[i].dividend * exp(-r * (q[i].time - t));
			}
		}
	}

	return SNon;
}

