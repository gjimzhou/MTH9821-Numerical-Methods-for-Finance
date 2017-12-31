//	Source file for Monte-Carlo Pricers

#include "MonteCarlo.hpp"

using namespace std;


//	Monte-Carlo Pricer for European Options.
vector<double> EuropeanMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(N, Generator, seed);

	double V = 0;
	double Delta = 0;
	double Vega = 0;

	//	Compute price and Greeks for each path.
	for (long int i = 0; i < N; i++)
	{
		double Si = S * exp((r - q - sig * sig / 2) * T + sig * sqrt(T) * z[i]);

		if (type == "C")
		{
			V += (exp(-r * T) * CallPayoff(Si, K));
			Delta += ((Si > K) * exp(-r * T) * Si / S);
			Vega += ((Si > K) * exp(-r * T) * Si * (-sig * T + sqrt(T) * z[i]));
		}
		else
		{
			V += (exp(-r * T) * PutPayoff(Si, K));
			Delta += (-(K > Si) * exp(-r * T) * Si / S);
			Vega += (-(K > Si) * exp(-r * T) * Si * (-sig * T + sqrt(T) * z[i]));
		}
	}

	//	Take average.
	V /= N;
	Delta /= N;
	Vega /= N;

	//	Output.
	vector<double> result;

	result.push_back(V);
	result.push_back(Delta);
	result.push_back(Vega);

	return result;
}

//	Monte-Carlo Pricer for Barrier Options.
double BarrierMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, double B, string Btype, long int n, long int m, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	long int N = n * m;
	vector<double> z = rnorm(N, Generator, seed);

	double dt = T / m;
	double V = 0;

	//	Generate n paths for underlying asset price.
	for (long int i = 0; i < n; i++)
	{
		//	Initial setup for option existence.
		bool in;
		if (Btype == "DnO" || Btype == "UnO") in = true;
		else in = false;

		//	Start with spot price.
		double Si = S;

		//	Generate m time intervals for each path.
		if (Btype == "UnI")
		{
			for (long int j = 0; j < m; j++)
			{
				Si = Si * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z[m * i + j]);

				//	If up then in.
				if (Si >= B) in = true;
			}
		}
		else if (Btype == "DnI")
		{
			for (long int j = 0; j < m; j++)
			{
				Si = Si * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z[m * i + j]);

				//	If down then in.
				if (Si <= B) in = true;
			}
		}
		else if (Btype == "UnO")
		{
			for (long int j = 0; j < m; j++)
			{
				Si = Si * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z[m * i + j]);

				//	If up then out.
				if (Si >= B) in = false;
			}
		}
		else 
		{
			for (long int j = 0; j < m; j++)
			{
				Si = Si * exp((r - q - sig * sig / 2) * dt + sig * sqrt(dt) * z[m * i + j]);

				//	If down then out.
				if (Si <= B) in = false;
			}
		}
		
		//	Compute price and Greeks for each path.
		if (type == "C")
		{
			V += (in * exp(-r * T) * CallPayoff(Si, K));
		}
		else
		{
			V += (in * exp(-r * T) * PutPayoff(Si, K));
		}
	}

	//	Take average.
	V /= n;

	return V;
}

//	Monte-Carlo Pricer for Barrier Options.	//	With optimal path number and time interval number.
double BarrierMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double q, double B, string Btype, long int N, string Generator, long int seed)
{
	long int m = ceil(pow(N, 1 / (double)3)*pow(T, 2 / (double)3));
	long int n = N / m;

	return BarrierMonteCarloPricer(type, S, K, T, r, sig, q, B, Btype, n, m, Generator, seed);
}


//	Monte-Carlo Pricer for European Options.	//	Control Variate Technique.
double EuropeanMCCVPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(N, Generator, seed);

	//	Initialize sample path.
	vector<double> ST;
	double Shat = 0;

	vector<double> V;
	double Vhat = 0;

	//	Compute stock price and payoff for each path.
	for (long int i = 0; i < N; i++)
	{
		double Si = S * exp((r - q - sig * sig / 2) * T + sig * sqrt(T) * z[i]);
		ST.push_back(Si);
		Shat += Si;

		if (type == "C")
		{
			double Vi = exp(-r * T) * CallPayoff(Si, K);
			V.push_back(Vi);
			Vhat += Vi;
		}
		else
		{
			double Vi = exp(-r * T) * PutPayoff(Si, K);
			V.push_back(Vi);
			Vhat += Vi;
		}
	}

	//	Take average.
	Shat /= N;
	Vhat /= N;

	//	Compute coefficient b.
	double Cov = 0;
	double Var = 0;

	for (long int i = 0; i < N; i++)
	{
		Cov += (ST[i] - Shat) * (V[i] - Vhat);
		Var += (ST[i] - Shat) * (ST[i] - Shat);
	}

	double b = Cov / Var;

	//	Compute payoff after control variate technique.
	double W = 0;

	for (long int i = 0; i < N; i++)
	{
		W += (V[i] - b * (ST[i] - exp(r * T) * S));
	}

	//	Take average.
	W /= N;

	return W;
}

//	Monte-Carlo Pricer for European Options.	//	Antithetic Variablews.
double EuropeanMCAVPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(N, Generator, seed);

	double V = 0;

	//	Compute stock price and payoff for each path.
	for (long int i = 0; i < N; i++)
	{
		double S1 = S * exp((r - q - sig * sig / 2) * T + sig * sqrt(T) * z[i]);
		double S2 = S * exp((r - q - sig * sig / 2) * T - sig * sqrt(T) * z[i]);

		if (type == "C")
		{
			V += (exp(-r * T) * CallPayoff(S1, K));
			V += (exp(-r * T) * CallPayoff(S2, K));
		}
		else
		{
			V += (exp(-r * T) * PutPayoff(S1, K));
			V += (exp(-r * T) * PutPayoff(S2, K));
		}
	}

	//	Take average.
	V /= (2 * N);

	return V;
}

//	Monte-Carlo Pricer for European Options.	//	Moment Matching.
double EuropeanMCMMPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(N, Generator, seed);

	//	Initialize sample path.
	vector<double> ST;
	double Shat = 0;
	
	double V = 0;

	//	Compute stock price for each path.
	for (long int i = 0; i < N; i++)
	{
		double Si = S * exp((r - q - sig * sig / 2) * T + sig * sqrt(T) * z[i]);
		ST.push_back(Si);
		Shat += Si;
	}

	//	Take average.
	Shat /= N;

	//	Compute payoff for each path.
	for (long int i = 0; i < N; i++)
	{
		//	Adjust stock price.
		ST[i] *= (exp(r * T) * S / Shat);

		if (type == "C")
		{
			V += exp(-r * T) * CallPayoff(ST[i], K);
		}
		else
		{
			V += exp(-r * T) * PutPayoff(ST[i], K);
		}
	}

	//	Take average.
	V /= N;

	return V;
}

//	Monte-Carlo Pricer for European Options.	//	Simultaneous Moment Matching and Control Variate Technique.
double EuropeanMCCVMMPricer(string type, double S, double K, double T, double r, double sig, double q, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(N, Generator, seed);

	//	Initialize sample path.
	vector<double> ST;
	double Shat = 0;

	vector<double> V;
	double Vhat = 0;

	//	Compute stock price for each path.
	for (long int i = 0; i < N; i++)
	{
		double Si = S * exp((r - q - sig * sig / 2) * T + sig * sqrt(T) * z[i]);
		ST.push_back(Si);
		Shat += Si;
	}

	//	Take average.
	Shat /= N;

	//	Compute payoff for each path.
	for (long int i = 0; i < N; i++)
	{
		//	Adjust stock price.
		ST[i] *= (exp(r * T) * S / Shat);

		if (type == "C")
		{
			double Vi = exp(-r * T) * CallPayoff(ST[i], K);
			V.push_back(Vi);
			Vhat += Vi;
		}
		else
		{
			double Vi = exp(-r * T) * PutPayoff(ST[i], K);
			V.push_back(Vi);
			Vhat += Vi;
		}
	}

	//	Take average.
	Shat = exp(r * T) * S;
	Vhat /= N;

	//	Compute coefficient b.
	double Cov = 0;
	double Var = 0;

	for (long int i = 0; i < N; i++)
	{
		Cov += (ST[i] - Shat) * (V[i] - Vhat);
		Var += (ST[i] - Shat) * (ST[i] - Shat);
	}

	double b = Cov / Var;

	//	Compute payoff after control variate technique.
	double W = 0;

	for (long int i = 0; i < N; i++)
	{
		W += (V[i] - b * (ST[i] - exp(r * T) * S));
	}

	//	Take average.
	W /= N;

	return W;
}


//	Monte-Carlo Pricer for European Basket Options.
double BasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int N, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	vector<double> z = rnorm(2 * N, Generator, seed);

	double V = 0;

	//	Compute price and Greeks for each path.
	for (long int i = 0; i < N; i++)
	{
		double S1i = S1 * exp((r - q1 - sig1 * sig1 / 2) * T + sig1 * sqrt(T) * z[2 * i]);
		double S2i = S2 * exp((r - q2 - sig2 * sig2 / 2) * T + sig2 * sqrt(T) * (rho * z[2 * i] + sqrt(1 - rho * rho) * z[2 * i + 1]));

		if (type == "C")
		{
			V += (exp(-r * T) * CallPayoff(S1i + S2i, K));
		}
		else
		{
			V += (exp(-r * T) * PutPayoff(S1i + S2i, K));
		}
	}

	//	Take average.
	V /= N;

	return V;
}

//	Monte-Carlo Pricer for Lookback Basket Options.
double LookbackBasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int n, long int m, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	long int N = 2 * n * m;
	vector<double> z = rnorm(N, Generator, seed);

	double dt = T / m;
	double V = 0;

	//	Generate n paths for underlying asset price.
	for (long int i = 0; i < n; i++)
	{
		//	Initial setup for maximum stock price.
		double Smax = S1 + S2;

		//	Start with spot price.
		double S1i = S1;
		double S2i = S2;

		//	Generate m time intervals for each path.
		for (long int j = 0; j < m; j++)
		{
			S1i = S1i * exp((r - q1 - sig1 * sig1 / 2) * dt + sig1 * sqrt(dt) * z[2 *(m * i + j)]);
			S2i = S2i * exp((r - q2 - sig2 * sig2 / 2) * dt + sig2 * sqrt(dt) * (rho * z[2 * (m * i + j)] + sqrt(1 - rho * rho) * z[2 * (m * i + j) + 1]));

			//	Update maximum stock price.
			if (S1i + S2i >= Smax) Smax = S1i + S2i;
		}

		//	Compute price and Greeks for each path.
		if (type == "C")
		{
			V += (exp(-r * T) * CallPayoff(Smax, K));
		}
		else
		{
			V += (exp(-r * T) * PutPayoff(Smax, K));
		}
	}

	//	Take average.
	V /= n;

	return V;
}

//	Monte-Carlo Pricer for Lookback Basket Options.	//	With optimal path number and time interval number.
double LookbackBasketMonteCarloPricer(string type, double S1, double S2, double K, double T, double r, double sig1, double sig2, double rho, double q1, double q2, long int N, string Generator, long int seed)
{
	long int m = ceil(pow(N, 1 / (double)3)*pow(T, 2 / (double)3));
	long int n = N / m;

	return LookbackBasketMonteCarloPricer(type, S1, S2, K, T, r, sig1, sig2, rho, q1, q2, n, m, Generator, seed);
}


//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.
double EuropeanHestonMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double lambda, double sigbar, double eta, double rho, double q, long int n, long int m, string Generator, long int seed)
{
	//	Generate standard normal distribution random variables.
	long int N = 2 * n * m;
	vector<double> z = rnorm(N, Generator, seed);

	double dt = T / m;
	double V = 0;

	//	Generate n paths for underlying asset price.
	for (long int i = 0; i < n; i++)
	{
		//	Start with spot price and initial variance.
		double Si = S;
		double sigi = sig;

		//	Generate m time intervals for each path.
		for (long int j = 0; j < m; j++)
		{
			Si = Si * exp((r - q - sigi * sigi / 2) * dt + sigi * sqrt(dt) * z[2 * (m * i + j)]);
			double sig2 = sigi * sigi - lambda * (sigi * sigi - sigbar * sigbar) * dt + eta * sigi * sqrt(dt) * (rho * z[2 * (m * i + j)] + sqrt(1 - rho * rho) * z[2 * (m * i + j) + 1]);
			sigi = sig2 > 0 ? sqrt(sig2) : 0;
		}

		//	Compute price and Greeks for each path.
		if (type == "C")
		{
			V += (exp(-r * T) * CallPayoff(Si, K));
		}
		else
		{
			V += (exp(-r * T) * PutPayoff(Si, K));
		}
	}

	//	Take average.
	V /= n;

	return V;
}

//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.	//	With optimal path number and time interval number.
double EuropeanHestonMonteCarloPricer(string type, double S, double K, double T, double r, double sig, double lambda, double sigbar, double eta, double rho, double q, long int N, string Generator, long int seed)
{
	long int m = ceil(pow(N, 1 / (double)3)*pow(T, 2 / (double)3));
	long int n = N / m;

	return EuropeanHestonMonteCarloPricer(type, S, K, T, r, sig, lambda, sigbar, eta, rho, q, n, m, Generator, seed);
}
