//	Source file for Binomial Pricers

#include "BinomialTree.hpp"

using namespace std;


//	Binomial Tree Pricer for European Options.
vector<double> EuropeanBinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (exp(-q * dt) - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - exp(-q * dt)) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	for (int k = N - 1; k >= 2; k--)
	{
		for (int i = 0; i <= k; i++)
		{
			V[i] = p1 * V[i] + p2 * V[i + 1];
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S * d;

	double S20 = S * u * u;
	double S21 = S;
	double S22 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10 = p1 * V20 + p2 * V21;
	double V11 = p1 * V21 + p2 * V22;

	double V00 = p1 * V10 + p2 * V11;

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}

//	Binomial Tree Pricer for American Options.
vector<double> AmericanBinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (exp(-q * dt) - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - exp(-q * dt)) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	if (type == "C")
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = S * pow(u, k - i) * pow(d, i);
				double payoff = CallPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}
	else
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = S * pow(u, k - i) * pow(d, i);
				double payoff = PutPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S * d;

	double S20 = S * u * u;
	double S21 = S;
	double S22 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10;
	double V11;
	double V00;

	if (type == "C")
	{
		double payoff10 = CallPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = CallPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = CallPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}
	else
	{
		double payoff10 = PutPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = PutPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = PutPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}


//	Average Binomial Tree Pricer for European Options.
vector<double> EuropeanABTPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	//	Binomial Tree Pricing for N and N+1.
	vector<double> VN1 = EuropeanBinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = EuropeanBinomialTreePricer(type, S, K, T, r, sig, q, N + 1);

	//	Average the prices and Greeks.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		result.push_back((VN1[i] + VN2[i]) / 2);
	}

	return result;
}

//	Average Binomial Tree Pricer for American Options.
vector<double> AmericanABTPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	//	Binomial Tree Pricing for N and N+1.
	vector<double> VN1 = AmericanBinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = AmericanBinomialTreePricer(type, S, K, T, r, sig, q, N + 1);

	//	Average the prices and Greeks.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		result.push_back((VN1[i] + VN2[i]) / 2);
	}

	return result;
}


//	Binomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanBBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (exp(-q * dt) - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - exp(-q * dt)) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity minus dt.
	if (type == "C")
	{
		for (int i = 0; i <= N - 1; i++)
		{
			double ST = S * pow(u, N - 1 - i) * pow(d, i);
			double VT = BSCallPrice(ST, K, dt, r, sig, q);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N - 1; i++)
		{
			double ST = S * pow(u, N - 1 - i) * pow(d, i);
			double VT = BSPutPrice(ST, K, dt, r, sig, q);
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	for (int k = N - 2; k >= 2; k--)
	{
		for (int i = 0; i <= k; i++)
		{
			V[i] = p1 * V[i] + p2 * V[i + 1];
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S * d;

	double S20 = S * u * u;
	double S21 = S;
	double S22 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10 = p1 * V20 + p2 * V21;
	double V11 = p1 * V21 + p2 * V22;

	double V00 = p1 * V10 + p2 * V11;

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}

//	Binomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanBBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (exp(-q * dt) - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - exp(-q * dt)) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity minus dt.
	if (type == "C")
	{
		for (int i = 0; i <= N - 1; i++)
		{
			double ST = S * pow(u, N - 1 - i) * pow(d, i);
			double priceT = BSCallPrice(ST, K, dt, r, sig, q);
			double payoffT = CallPayoff(ST, K);
			double VT = priceT > payoffT ? priceT : payoffT;
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N - 1; i++)
		{
			double ST = S * pow(u, N - 1 - i) * pow(d, i);
			double priceT = BSPutPrice(ST, K, dt, r, sig, q);
			double payoffT = PutPayoff(ST, K);
			double VT = priceT > payoffT ? priceT : payoffT;
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	if (type == "C")
	{
		for (int k = N - 2; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = S * pow(u, k - i) * pow(d, i);
				double payoff = CallPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}
	else
	{
		for (int k = N - 2; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = S * pow(u, k - i) * pow(d, i);
				double payoff = PutPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S * d;

	double S20 = S * u * u;
	double S21 = S;
	double S22 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10;
	double V11;
	double V00;

	if (type == "C")
	{
		double payoff10 = CallPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = CallPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = CallPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}
	else
	{
		double payoff10 = PutPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = PutPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = PutPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}


//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanBBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	int M = N / 2;

	//	Binomial Tree Black-Scholes Pricing for N and N/2.
	vector<double> VN1 = EuropeanBBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = EuropeanBBSPricer(type, S, K, T, r, sig, q, M);

	//	Caculate the prices and Greeks with Richardson Extrapolation.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		double res = 2 * VN1[i] - VN2[i];
		result.push_back(res);
	}

	return result;
}

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanBBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	int M = N / 2;

	//	Binomial Tree Black-Scholes Pricing for N and N+1.
	vector<double> VN1 = AmericanBBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = AmericanBBSPricer(type, S, K, T, r, sig, q, M);

	//	Caculate the prices and Greeks with Richardson Extrapolation.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		double res = 2 * VN1[i] - VN2[i];
		result.push_back(res);
	}

	return result;
}


//	Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanBinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanBinomialTreePricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Average Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanABTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanABTPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanABTPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Binomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBBSVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanBBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanBBSPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanBBSRVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanBBSRPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanBBSRPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}


//	Implied Volatility for American Options using Binomial Tree Pricer.
double AmericanBTImpVol(string type, double Vm, double S, double K, double T, double r, double q, int N, double sig0, double sig1, double thres)
{
	while (abs(sig1 - sig0) > thres)
	{
		vector<double> pricer0 = AmericanBinomialTreePricer(type, S, K, T, r, sig0, q, N);
		double V0 = pricer0[0];
		vector<double> pricer1 = AmericanBinomialTreePricer(type, S, K, T, r, sig1, q, N);
		double V1 = pricer1[0];

		double sigtmp = sig1;
		sig1 = sig1 - (V1 - Vm) * (sig1 - sig0) / (V1 - V0);
		sig0 = sigtmp;

		cout << sig1 << endl;
	}

	return sig1;
}




//	Binomial Tree Pricer for European Options.
vector<double> EuropeanDividendBinomialTreePricer(string type, double S, double K, double T, double r, double sig, vector<Dividend> q, int N)
{
	S = DivToNon(q, S, 0, T, r);
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (1 - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - 1) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	for (int k = N - 1; k >= 2; k--)
	{
		for (int i = 0; i <= k; i++)
		{
			V[i] = p1 * V[i] + p2 * V[i + 1];
		}
	}

	//	Final two steps for pricing.
	double S10 = NonToDiv(q, S * u, 1 * dt, T, r);
	double S11 = NonToDiv(q, S * d, 1 * dt, T, r);

	double S20 = NonToDiv(q, S * u * u, 1 * dt, T, r);
	double S21 = NonToDiv(q, S, 1 * dt, T, r);
	double S22 = NonToDiv(q, S * d * d, 1 * dt, T, r);

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10 = p1 * V20 + p2 * V21;
	double V11 = p1 * V21 + p2 * V22;

	double V00 = p1 * V10 + p2 * V11;

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}

//	Binomial Tree Pricer for American Options.
vector<double> AmericanDividendBinomialTreePricer(string type, double S, double K, double T, double r, double sig, vector<Dividend> q, int N)
{
	S = DivToNon(q, S, 0, T, r);
	double dt = T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1 / u;

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = (1 - d * exp(-r * dt)) / (u - d);
	double p2 = (u * exp(-r * dt) - 1) / (u - d);

	//	Equity binomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= N; i++)
		{
			double ST = S * pow(u, N - i) * pow(d, i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Binomial Tree Pricing until T=2.
	if (type == "C")
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = NonToDiv(q, S * pow(u, k - i) * pow(d, i), k * dt, T, r);
				double payoff = CallPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}
	else
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1];
				double St = NonToDiv(q, S * pow(u, k - i) * pow(d, i), k * dt, T, r);
				double payoff = PutPayoff(St, K);
				//	if (payoff > price) cout << "Early Ex!" << endl;
				V[i] = price > payoff ? price : payoff;
			}
		}
	}

	//	Final two steps for pricing.
	double S10 = NonToDiv(q, S * u, 1 * dt, T, r);
	double S11 = NonToDiv(q, S * d, 1 * dt, T, r);

	double S20 = NonToDiv(q, S * u * u, 1 * dt, T, r);
	double S21 = NonToDiv(q, S, 1 * dt, T, r);
	double S22 = NonToDiv(q, S * d * d, 1 * dt, T, r);

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];

	double V10;
	double V11;
	double V00;

	if (type == "C")
	{
		double payoff10 = CallPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = CallPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = CallPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}
	else
	{
		double payoff10 = PutPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = PutPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff00 = PutPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}

	//	Calculate Greeks.
	double Delta00 = (V10 - V11) / (S10 - S11);
	double Delta10 = (V20 - V21) / (S20 - S21);
	double Delta11 = (V21 - V22) / (S21 - S22);

	double Gamma00 = 2 * (Delta10 - Delta11) / (S20 - S22);

	double Theta00 = (V21 - V00) / (2 * dt);

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}
