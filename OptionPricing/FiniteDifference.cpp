//	Source file for Finite Difference Pricers

#include "FiniteDifference.hpp"

using namespace std;


//	Finite Difference Pricer for European Options.
vector<double> EuropeanExplicitFDPricer(string type, double S, double K, double T, double r, double sig, double q, long int M, double alpha)
{
	//	Black-Scholes price.
	double VBS;
	if (type == "C") VBS = BSCallPrice(S, K, T, r, sig, q);
	else VBS = BSPutPrice(S, K, T, r, sig, q);

	//	Change of variables constants.
	double a = (r - q) / (sig * sig) - 0.5;
	double b = ((r - q) / (sig * sig) + 0.5) * ((r - q) / (sig * sig) + 0.5) + 2 * q / (sig * sig);

	//	Computational domain.
	double taufinal = T * sig * sig / 2;
	double xleft = log(S / K) + (r - q - sig * sig / 2) * T - 3 * sig * sqrt(T);
	double xright = log(S / K) + (r - q - sig * sig / 2) * T + 3 * sig * sqrt(T);

	//	Discretization.
	double dtau = taufinal / M;
	double dx = sqrt(dtau / alpha);
	long int N = int(floor((xright - xleft) / dx));
	dx = (xright - xleft) / N;
	alpha = dtau / (dx * dx);

	//	Build scheme.
	vector<double> tau;
	vector<double> x;

	for (long int m = 0; m < M + 1; m++)
	{
		tau.push_back(m * dtau);
	}

	for (long int n = 0; n < N + 1; n++)
	{
		x.push_back(xleft + n * dx);
	}

	//	Heat equation matrix.
	vector<vector<double>> u;
	for (long int m = 0; m < M + 1; m++)
	{
		u.push_back(x);
	}

	//	Initial condition.
	for (long int n = 0; n < N + 1; n++)
	{
		double upayoff;
		if (type == "C") upayoff = (exp(x[n]) - 1) > 0 ? (exp(x[n]) - 1) : 0;
		else upayoff = (1 - exp(x[n])) > 0 ? (1 - exp(x[n])) : 0;

		u[0][n] = K * exp(a * x[n]) * upayoff;
	}

	//	Boundary condition.
	for (long int m = 1; m < M + 1; m++)
	{
		if (type == "C")
		{
			u[m][0] = 0;
			u[m][N] = K * exp(a * xright + b * tau[m]) * (exp(xright - 2 * q * tau[m] / (sig * sig)) - exp(-2 * r * tau[m] / (sig * sig)));
		}
		else
		{
			u[m][0] = K * exp(a * xleft + b * tau[m]) * (exp(-2 * r * tau[m] / (sig * sig)) - exp(xleft - 2 * q * tau[m] / (sig * sig)));
			u[m][N] = 0;
		}
	}

	//	Finite Differences.
	for (long int m = 0; m < M; m++)
	{
		for (long int n = 1; n < N; n++)
		{
			u[m + 1][n] = alpha * u[m][n - 1] + (1 - 2 * alpha) * u[m][n] + alpha * u[m][n + 1];
		}
	}

	//	Price at each point.
	vector<double> S0;
	for (long int n = 0; n < N + 1; n++)
	{
		S0.push_back(K * exp(x[n]));
	}

	vector<double> Vex;
	for (long int n = 0; n < N + 1; n++)
	{
		if (type == "C") Vex.push_back(BSCallPrice(S0[n], K, T, r, sig, q));
		else Vex.push_back(BSPutPrice(S0[n], K, T, r, sig, q));
	}

	vector<double> Vapp;
	for (long int n = 0; n < N + 1; n++)
	{
		Vapp.push_back(exp(-a * x[n] - b * taufinal) * u[M][n]);
	}

	//	Pointwise value.
	double xcom = log(S / K);

	long int i = 0;
	while (x[i] < xcom)
	{
		i++;
	}
	i--;

	double Vapp1 = ((S0[i + 1] - S) * Vapp[i] + (S - S0[i]) * Vapp[i + 1]) / (S0[i + 1] - S0[i]);
	double Vapp2 = exp(-a * xcom - b * taufinal) * ((x[i + 1] - xcom) * u[M][i] + (xcom - x[i]) * u[M][i + 1]) / (x[i + 1] - x[i]);

	//	Pointwise error.
	double errpt1 = abs(Vapp1 - VBS);
	double errpt2 = abs(Vapp2 - VBS);

	//	RMS error.
	double errRMS = 0;
	long int NRMS = 0;

	for (long int n = 0; n < N + 1; n++)
	{
		if (Vex[n] > 0.00001 * S0[n])
		{
			errRMS += (((Vapp[n] - Vex[n]) * (Vapp[n] - Vex[n])) / (Vex[n] * Vex[n]));
			NRMS++;
		}
	}

	errRMS /= NRMS;

	//	Greeks.
	double Sim1 = K * exp(x[i - 1]);
	double Si2 = K * exp(x[i + 2]);

	double Delta = (Vapp[i + 1] - Vapp[i]) / (S0[i + 1] - S0[i]);
	double Gamma = ((Vapp[i + 2] - Vapp[i + 1]) / (S0[i + 2] - S0[i + 1]) - (Vapp[i] - Vapp[i - 1]) / (S0[i] - S0[i - 1])) / ((S0[i + 2] + S0[i + 1]) / 2 - (S0[i] + S0[i - 1]) / 2);

	double dt = 2 * dtau / (sig * sig);
	double Vi0dt = exp(-a * x[i] - b * (taufinal - dtau)) * u[M - 1][i];
	double Vi1dt = exp(-a * x[i + 1] - b * (taufinal - dtau)) * u[M - 1][i + 1];
	double Vapp1dt = ((S0[i + 1] - S) * Vi0dt + (S - S0[i]) * Vi1dt) / (S0[i + 1] - S0[i]);

	double Theta = (Vapp1dt - Vapp1) / dt;

	//	Output.
	vector<double> res;
	res.push_back(Vapp1);
	res.push_back(errpt1);
	res.push_back(Vapp2);
	res.push_back(errpt2);
	res.push_back(errRMS);
	res.push_back(Delta);
	res.push_back(Gamma);
	res.push_back(Theta);	

	return res;
}
