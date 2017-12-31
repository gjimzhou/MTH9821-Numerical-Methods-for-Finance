//	Source file for Homework 1

#include <iostream>

#include "HW1.hpp"

using namespace std;


//	Exercise 1: Derivatives Valuation and Hedging Parameters Computation.

void HW1_Test1E(double& S, EuropeanOption& EOSample, ofstream& file)
{
	file << "Binomial Tree Methods for European Options" << endl;

	for (int N = 10; N <= 100; N++)
	{
		vector<double> BT = EOSample.BinomialTreePricer(S, N);
		vector<double> ABT = EOSample.ABTPricer(S, N);
		vector<double> BBS = EOSample.BBSPricer(S, N);
		vector<double> BBSR = EOSample.BBSRPricer(S, N);

		file << BT[0] << "," << ABT[0] << "," << BBS[0] << "," << BBSR[0] << endl;
	}

	file << endl;
}

void HW1_Test1A(double& S, AmericanOption& AOSample, ofstream& file)
{
	file << "Binomial Tree Methods for American Options" << endl;

	for (int N = 10; N <= 100; N++)
	{
		vector<double> BT = AOSample.BinomialTreePricer(S, N);
		vector<double> ABT = AOSample.ABTPricer(S, N);
		vector<double> BBS = AOSample.BBSPricer(S, N);
		vector<double> BBSR = AOSample.BBSRPricer(S, N);

		file << BT[0] << "," << ABT[0] << "," << BBS[0] << "," << BBSR[0] << endl;
	}

	file << endl;
}

void HW1_Test1(double& S, EuropeanOption& EOSample, AmericanOption& AOSample, ofstream& file)
{
	file << "Exercise 1: Binomial Tree Methods for Derivatives Valuation and Hedging Parameters Computation." << endl;

	HW1_Test1E(S, EOSample, file);
	HW1_Test1A(S, AOSample, file);

	file << endl;
}


//	Exercise 2: European Option.

void HW1_Test2BT(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Binomial Tree for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BT = EOSample.BinomialTreePricer(S, N);

		double V1 = BT[0];
		double Delta1 = BT[1];
		double Gamma1 = BT[2];
		double Theta1 = BT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test2ABT(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Average Binomial Tree for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> ABT = EOSample.ABTPricer(S, N);

		double V1 = ABT[0];
		double Delta1 = ABT[1];
		double Gamma1 = ABT[2];
		double Theta1 = ABT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test2BBS(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Binomial Tree Black-Scholes for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBS = EOSample.BBSPricer(S, N);

		double V1 = BBS[0];
		double Delta1 = BBS[1];
		double Gamma1 = BBS[2];
		double Theta1 = BBS[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test2BBSR(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Binomial Tree Black-Scholes with Richardson Extrapolation for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBSR = EOSample.BBSRPricer(S, N);

		double V1 = BBSR[0];
		double Delta1 = BBSR[1];
		double Gamma1 = BBSR[2];
		double Theta1 = BBSR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test2(double& S, EuropeanOption& EOSample, ofstream& file)
{
	file << "Exercise 2: Binomial Tree Methods for European Options." << endl;

	vector<double> BS = EOSample.BSPricer(S);

	HW1_Test2BT(S, EOSample, BS, file);
	HW1_Test2ABT(S, EOSample, BS, file);
	HW1_Test2BBS(S, EOSample, BS, file);
	HW1_Test2BBSR(S, EOSample, BS, file);

	file << endl;
}


//	Exercise 3: American Option.

void HW1_Test3BT(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BT = AOSample.BinomialTreePricer(S, N);

		double V1 = BT[0];
		double Delta1 = BT[1];
		double Gamma1 = BT[2];
		double Theta1 = BT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3ABT(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Average Binomial Tree for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> ABT = AOSample.ABTPricer(S, N);

		double V1 = ABT[0];
		double Delta1 = ABT[1];
		double Gamma1 = ABT[2];
		double Theta1 = ABT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3BBS(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree Black-Scholes for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBS = AOSample.BBSPricer(S, N);

		double V1 = BBS[0];
		double Delta1 = BBS[1];
		double Gamma1 = BBS[2];
		double Theta1 = BBS[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3BBSR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree Black-Scholes with Richardson Extrapolation for American Options" << endl;
	
	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBSR = AOSample.BBSRPricer(S, N);

		double V1 = BBSR[0];
		double Delta1 = BBSR[1];
		double Gamma1 = BBSR[2];
		double Theta1 = BBSR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3BTVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BTVR = AOSample.BTVRPricer(S, N);

		double V1 = BTVR[0];
		double Delta1 = BTVR[1];
		double Gamma1 = BTVR[2];
		double Theta1 = BTVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3ABTVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Average Binomial Tree with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> ABTVR = AOSample.ABTVRPricer(S, N);

		double V1 = ABTVR[0];
		double Delta1 = ABTVR[1];
		double Gamma1 = ABTVR[2];
		double Theta1 = ABTVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3BBSVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree Black-Scholes with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBSVR = AOSample.BBSVRPricer(S, N);

		double V1 = BBSVR[0];
		double Delta1 = BBSVR[1];
		double Gamma1 = BBSVR[2];
		double Theta1 = BBSVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3BBSRVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Binomial Tree Black-Scholes with Richardson Extrapolation with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> BBSRVR = AOSample.BBSRVRPricer(S, N);

		double V1 = BBSRVR[0];
		double Delta1 = BBSRVR[1];
		double Gamma1 = BBSRVR[2];
		double Theta1 = BBSRVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW1_Test3(double& S, AmericanOption& AOSample, ofstream& file)
{
	file << "Exercise 3: Binomial Tree Methods for American Options." << endl;

	vector<double> Exact = AOSample.ABTPricer(S, 10000);

	HW1_Test3BT(S, AOSample, Exact, file);
	HW1_Test3ABT(S, AOSample, Exact, file);
	HW1_Test3BBS(S, AOSample, Exact, file);
	HW1_Test3BBSR(S, AOSample, Exact, file);

	HW1_Test3BTVR(S, AOSample, Exact, file);
	HW1_Test3ABTVR(S, AOSample, Exact, file);
	HW1_Test3BBSVR(S, AOSample, Exact, file);
	HW1_Test3BBSRVR(S, AOSample, Exact, file);

	file << endl;
}


void HW1_Test()
{
	ofstream file("HW1.csv");
	file << "!! MTH9821 Homework 1 !!" << endl;

	int precision = numeric_limits<double>::max_digits10;
	file << setprecision(precision) << endl;

	double S = 45.0;
	EuropeanOption EOSample(41.0, 1.0, 0.025, 0.35, 0.01, "P", "EOSample");
	AmericanOption AOSample(41.0, 1.0, 0.025, 0.35, 0.01, "P", "AOSample");

	HW1_Test1(S, EOSample, AOSample, file);
	HW1_Test2(S, EOSample, file);
	HW1_Test3(S, AOSample, file);

	file << endl;
}

