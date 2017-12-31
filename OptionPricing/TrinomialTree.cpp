//	Source file for Trinomial Pricers

#include "TrinomialTree.hpp"

using namespace std;


//	Trinomial Tree Pricer for European Options.
vector<double> EuropeanTrinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1 / u;
	double dis = exp(-r * dt);

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = dis * (1 / (double)6 + (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));
	double p2 = dis * 2 / 3;
	double p3 = dis * (1 / (double)6 - (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));

	//	Equity Trinomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= 2 * N; i++)
		{
			double ST = S * pow(u, N - i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= 2 * N; i++)
		{
			double ST = S * pow(u, N - i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Trinomial Tree Pricing until T=2.
	for (int k = N - 1; k >= 2; k--)
	{
		for (int i = 0; i <= 2 * k; i++)
		{
			V[i] = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S;
	double S12 = S * d;

	double S20 = S * u * u;
	double S21 = S * u;
	double S22 = S;
	double S23 = S * d;
	double S24 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];
	double V23 = V[3];
	double V24 = V[4];

	double V10 = p1 * V20 + p2 * V21 + p3 * V22;
	double V11 = p1 * V21 + p2 * V22 + p3 * V23;
	double V12 = p1 * V22 + p2 * V23 + p3 * V24;

	double V00 = p1 * V10 + p2 * V11 + p3 * V12;

	//	Calculate Greeks.
	double Delta00 = (V10 - V12) / (S10 - S12);
	double Delta10 = (V20 - V22) / (S20 - S22);
	double Delta12 = (V22 - V24) / (S22 - S24);

	double Gamma00 = (Delta10 - Delta12) / (S10 - S12);

	double Theta00 = (V11 - V00) / dt;

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}

//	Trinomial Tree Pricer for American Options.
vector<double> AmericanTrinomialTreePricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1 / u;
	double dis = exp(-r * dt);

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = dis * (1 / (double)6 + (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));
	double p2 = dis * 2 / 3;
	double p3 = dis * (1 / (double)6 - (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));

	//	Equity Trinomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= 2 * N; i++)
		{
			double ST = S * pow(u, N - i);
			double VT = CallPayoff(ST, K);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= 2 * N; i++)
		{
			double ST = S * pow(u, N - i);
			double VT = PutPayoff(ST, K);
			V.push_back(VT);
		}
	}

	//	Backward Trinomial Tree Pricing until T=2.
	if (type == "C")
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= 2 * k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
				double St = S * pow(u, k - i);
				double payoff = CallPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}
	else
	{
		for (int k = N - 1; k >= 2; k--)
		{
			for (int i = 0; i <= 2 * k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
				double St = S * pow(u, k - i);
				double payoff = PutPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S;
	double S12 = S * d;

	double S20 = S * u * u;
	double S21 = S * u;
	double S22 = S;
	double S23 = S * d;
	double S24 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];
	double V23 = V[3];
	double V24 = V[4];

	double V10;
	double V11;
	double V12;
	double V00;

	if (type == "C")
	{
		double payoff10 = CallPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21 + p3 * V22;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = CallPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22 + p3 * V23;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff12 = CallPayoff(S12, K);
		double price12 = p1 * V22 + p2 * V23 + p3 * V24;
		V12 = price12 > payoff12 ? price12 : payoff12;

		double payoff00 = CallPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11 + p3 * V12;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}
	else
	{
		double payoff10 = PutPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21 + p3 * V22;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = PutPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22 + p3 * V23;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff12 = PutPayoff(S12, K);
		double price12 = p1 * V22 + p2 * V23 + p3 * V24;
		V12 = price12 > payoff12 ? price12 : payoff12;

		double payoff00 = PutPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11 + p3 * V12;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}

	//	Calculate Greeks.
	double Delta00 = (V10 - V12) / (S10 - S12);
	double Delta10 = (V20 - V22) / (S20 - S22);
	double Delta12 = (V22 - V24) / (S22 - S24);

	double Gamma00 = (Delta10 - Delta12) / (S10 - S12);

	double Theta00 = (V11 - V00) / dt;

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}


//	Average Trinomial Tree Pricer for European Options.
vector<double> EuropeanATTPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	//	Trinomial Tree Pricing for N and N+1.
	vector<double> VN1 = EuropeanTrinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = EuropeanTrinomialTreePricer(type, S, K, T, r, sig, q, N + 1);

	//	Average the prices and Greeks.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		result.push_back((VN1[i] + VN2[i]) / 2);
	}

	return result;
}

//	Average Trinomial Tree Pricer for American Options.
vector<double> AmericanATTPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	//	Trinomial Tree Pricing for N and N+1.
	vector<double> VN1 = AmericanTrinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = AmericanTrinomialTreePricer(type, S, K, T, r, sig, q, N + 1);

	//	Average the prices and Greeks.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		result.push_back((VN1[i] + VN2[i]) / 2);
	}

	return result;
}


//	Trinomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanTBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1 / u;
	double dis = exp(-r * dt);

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = dis * (1 / (double)6 + (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));
	double p2 = dis * 2 / 3;
	double p3 = dis * (1 / (double)6 - (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));

	//	Equity Trinomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= 2 * (N - 1); i++)
		{
			double ST = S * pow(u, N - 1 - i);
			double VT = BSCallPrice(ST, K, dt, r, sig, q);
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= 2 * (N - 1); i++)
		{
			double ST = S * pow(u, N - 1 - i);
			double VT = BSPutPrice(ST, K, dt, r, sig, q);
			V.push_back(VT);
		}
	}

	//	Backward Trinomial Tree Pricing until T=2.
	for (int k = N - 2; k >= 2; k--)
	{
		for (int i = 0; i <= 2 * k; i++)
		{
			V[i] = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S;
	double S12 = S * d;

	double S20 = S * u * u;
	double S21 = S * u;
	double S22 = S;
	double S23 = S * d;
	double S24 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];
	double V23 = V[3];
	double V24 = V[4];

	double V10 = p1 * V20 + p2 * V21 + p3 * V22;
	double V11 = p1 * V21 + p2 * V22 + p3 * V23;
	double V12 = p1 * V22 + p2 * V23 + p3 * V24;

	double V00 = p1 * V10 + p2 * V11 + p3 * V12;

	//	Calculate Greeks.
	double Delta00 = (V10 - V12) / (S10 - S12);
	double Delta10 = (V20 - V22) / (S20 - S22);
	double Delta12 = (V22 - V24) / (S22 - S24);

	double Gamma00 = (Delta10 - Delta12) / (S10 - S12);

	double Theta00 = (V11 - V00) / dt;

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}

//	Trinomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanTBSPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	double dt = T / N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1 / u;
	double dis = exp(-r * dt);

	//	Risk-neutral probability with discounted to simplify computation.
	double p1 = dis * (1 / (double)6 + (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));
	double p2 = dis * 2 / 3;
	double p3 = dis * (1 / (double)6 - (r - q - sig * sig / 2) * sqrt(dt / (12 * sig * sig)));

	//	Equity Trinomial tree value.
	vector<double> V;

	// Value at maturity.
	if (type == "C")
	{
		for (int i = 0; i <= 2 * (N - 1); i++)
		{
			double ST = S * pow(u, N - 1 - i);
			double priceT = BSCallPrice(ST, K, dt, r, sig, q);
			double payoffT = CallPayoff(ST, K);
			double VT = priceT > payoffT ? priceT : payoffT;
			V.push_back(VT);
		}
	}
	else
	{
		for (int i = 0; i <= 2 * (N - 1); i++)
		{
			double ST = S * pow(u, N - 1 - i);
			double priceT = BSPutPrice(ST, K, dt, r, sig, q);
			double payoffT = PutPayoff(ST, K);
			double VT = priceT > payoffT ? priceT : payoffT;
			V.push_back(VT);
		}
	}

	//	Backward Trinomial Tree Pricing until T=2.
	if (type == "C") 
	{
		for (int k = N - 2; k >= 2; k--)
		{
			for (int i = 0; i <= 2 * k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
				double St = S * pow(u, k - i);
				double payoff = CallPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}
	else
	{
		for (int k = N - 2; k >= 2; k--)
		{
			for (int i = 0; i <= 2 * k; i++)
			{
				double price = p1 * V[i] + p2 * V[i + 1] + p3 * V[i + 2];
				double St = S * pow(u, k - i);
				double payoff = PutPayoff(St, K);
				V[i] = price > payoff ? price : payoff;
			}
		}
	}

	//	Final two steps for pricing.
	double S10 = S * u;
	double S11 = S;
	double S12 = S * d;

	double S20 = S * u * u;
	double S21 = S * u;
	double S22 = S;
	double S23 = S * d;
	double S24 = S * d * d;

	double V20 = V[0];
	double V21 = V[1];
	double V22 = V[2];
	double V23 = V[3];
	double V24 = V[4];

	double V10;
	double V11;
	double V12;
	double V00;

	if (type == "C")
	{
		double payoff10 = CallPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21 + p3 * V22;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = CallPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22 + p3 * V23;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff12 = CallPayoff(S12, K);
		double price12 = p1 * V22 + p2 * V23 + p3 * V24;
		V12 = price12 > payoff12 ? price12 : payoff12;

		double payoff00 = CallPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11 + p3 * V12;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}
	else
	{
		double payoff10 = PutPayoff(S10, K);
		double price10 = p1 * V20 + p2 * V21 + p3 * V22;
		V10 = price10 > payoff10 ? price10 : payoff10;

		double payoff11 = PutPayoff(S11, K);
		double price11 = p1 * V21 + p2 * V22 + p3 * V23;
		V11 = price11 > payoff11 ? price11 : payoff11;

		double payoff12 = PutPayoff(S12, K);
		double price12 = p1 * V22 + p2 * V23 + p3 * V24;
		V12 = price12 > payoff12 ? price12 : payoff12;

		double payoff00 = PutPayoff(S, K);
		double price00 = p1 * V10 + p2 * V11 + p3 * V12;
		V00 = price00 > payoff00 ? price00 : payoff00;
	}

	//	Calculate Greeks.
	double Delta00 = (V10 - V12) / (S10 - S12);
	double Delta10 = (V20 - V22) / (S20 - S22);
	double Delta12 = (V22 - V24) / (S22 - S24);

	double Gamma00 = (Delta10 - Delta12) / (S10 - S12);

	double Theta00 = (V11 - V00) / dt;

	//	Output.
	vector<double> result;

	result.push_back(V00);
	result.push_back(Delta00);
	result.push_back(Gamma00);
	result.push_back(Theta00);

	return result;
}


//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanTBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	int M = N / 2;

	//	Trinomial Tree Black-Scholes Pricing for N and N/2.
	vector<double> VN1 = EuropeanTBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = EuropeanTBSPricer(type, S, K, T, r, sig, q, M);

	//	Caculate the prices and Greeks with Richardson Extrapolation.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		double res = 2 * VN1[i] - VN2[i];
		result.push_back(res);
	}

	return result;
}

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanTBSRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	int M = N / 2;

	//	Trinomial Tree Black-Scholes Pricing for N and N+1.
	vector<double> VN1 = AmericanTBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> VN2 = AmericanTBSPricer(type, S, K, T, r, sig, q, M);

	//	Caculate the prices and Greeks with Richardson Extrapolation.
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		double res = 2 * VN1[i] - VN2[i];
		result.push_back(res);
	}

	return result;
}


//	Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanTrinomialTreePricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanTrinomialTreePricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Average Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanATTVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanATTPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanATTPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Trinomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTBSVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanTBSPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanTBSPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanTBSRVRPricer(string type, double S, double K, double T, double r, double sig, double q, int N)
{
	vector<double> BS = EuropeanBSPricer(type, S, K, T, r, sig, q);
	vector<double> EO = EuropeanTBSRPricer(type, S, K, T, r, sig, q, N);
	vector<double> AO = AmericanTBSRPricer(type, S, K, T, r, sig, q, N);

	vector<double> result;

	//	Variance reducgtion by subtract the error from European to BS.
	for (int i = 0; i < 4; i++)
	{
		result.push_back(AO[i] - (EO[i] - BS[i]));
	}

	return result;
}

