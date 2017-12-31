//	Source file for Homework 2

#include <iostream>

#include "HW2.hpp"

using namespace std;


//	Exercise 1: European Option.

void HW2_Test1TT(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Trinomial Tree for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TT = EOSample.TrinomialTreePricer(S, N);

		double V1 = TT[0];
		double Delta1 = TT[1];
		double Gamma1 = TT[2];
		double Theta1 = TT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test1ATT(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Average Trinomial Tree for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> ATT = EOSample.ATTPricer(S, N);

		double V1 = ATT[0];
		double Delta1 = ATT[1];
		double Gamma1 = ATT[2];
		double Theta1 = ATT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test1TBS(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBS = EOSample.TBSPricer(S, N);

		double V1 = TBS[0];
		double Delta1 = TBS[1];
		double Gamma1 = TBS[2];
		double Theta1 = TBS[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test1TBSR(double& S, EuropeanOption& EOSample, vector<double>& BS, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes with Richardson Extrapolation for European Options" << endl;

	double V = BS[0];
	double Delta = BS[1];
	double Gamma = BS[2];
	double Theta = BS[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBSR = EOSample.TBSRPricer(S, N);

		double V1 = TBSR[0];
		double Delta1 = TBSR[1];
		double Gamma1 = TBSR[2];
		double Theta1 = TBSR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test1(double& S, EuropeanOption& EOSample, ofstream& file)
{
	file << "Exercise 1: Trinomial Tree Methods for European Options." << endl;

	vector<double> BS = EOSample.BSPricer(S);

	HW2_Test1TT(S, EOSample, BS, file);
	HW2_Test1ATT(S, EOSample, BS, file);
	HW2_Test1TBS(S, EOSample, BS, file);
	HW2_Test1TBSR(S, EOSample, BS, file);

	file << endl;
}


//	Exercise 2: American Option.

void HW2_Test2TT(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TT = AOSample.TrinomialTreePricer(S, N);

		double V1 = TT[0];
		double Delta1 = TT[1];
		double Gamma1 = TT[2];
		double Theta1 = TT[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2TBS(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBS = AOSample.TBSPricer(S, N);

		double V1 = TBS[0];
		double Delta1 = TBS[1];
		double Gamma1 = TBS[2];
		double Theta1 = TBS[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2TBSR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes with Richardson Extrapolation for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBSR = AOSample.TBSRPricer(S, N);

		double V1 = TBSR[0];
		double Delta1 = TBSR[1];
		double Gamma1 = TBSR[2];
		double Theta1 = TBSR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2TTVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TTVR = AOSample.TTVRPricer(S, N);

		double V1 = TTVR[0];
		double Delta1 = TTVR[1];
		double Gamma1 = TTVR[2];
		double Theta1 = TTVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2TBSVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBSVR = AOSample.TBSVRPricer(S, N);

		double V1 = TBSVR[0];
		double Delta1 = TBSVR[1];
		double Gamma1 = TBSVR[2];
		double Theta1 = TBSVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2TBSRVR(double& S, AmericanOption& AOSample, vector<double>& Exact, ofstream& file)
{
	file << "Trinomial Tree Black-Scholes with Richardson Extrapolation with variance reduction for American Options" << endl;

	double V = Exact[0];
	double Delta = Exact[1];
	double Gamma = Exact[2];
	double Theta = Exact[3];

	file << V << "," << Delta << "," << Gamma << "," << Theta << endl;

	for (int N = 10; N <= 1280; N *= 2)
	{
		vector<double> TBSRVR = AOSample.TBSRVRPricer(S, N);

		double V1 = TBSRVR[0];
		double Delta1 = TBSRVR[1];
		double Gamma1 = TBSRVR[2];
		double Theta1 = TBSRVR[3];

		double VError = abs(V1 - V);
		double DeltaError = abs(Delta1 - Delta);
		double GammaError = abs(Gamma1 - Gamma);
		double ThetaError = abs(Theta1 - Theta);

		file << N << "," << V1 << "," << VError << "," << N * VError << "," << N * N * VError << "," <<
			Delta1 << "," << DeltaError << "," << Gamma1 << "," << GammaError << "," << Theta1 << "," << ThetaError << endl;
	}

	file << endl;
}

void HW2_Test2(double& S, AmericanOption& AOSample, ofstream& file)
{
	file << "Exercise 2: Trinomial Tree Methods for American Options." << endl;

	vector<double> Exact = AOSample.ABTPricer(S, 10000);

	HW2_Test2TT(S, AOSample, Exact, file);
	HW2_Test2TBS(S, AOSample, Exact, file);
	HW2_Test2TBSR(S, AOSample, Exact, file);

	HW2_Test2TTVR(S, AOSample, Exact, file);
	HW2_Test2TBSVR(S, AOSample, Exact, file);
	HW2_Test2TBSRVR(S, AOSample, Exact, file);

	file << endl;
}


void HW2_Test()
{
	ofstream file("HW2.csv");
	file << "!! MTH9821 Homework 2 !!" << endl;

	int precision = numeric_limits<double>::max_digits10;
	file << setprecision(precision) << endl;

	double S = 41.0;
	EuropeanOption EOSample(39.0, 1.0, 0.03, 0.25, 0.005, "P", "EOSample");
	AmericanOption AOSample(39.0, 1.0, 0.03, 0.25, 0.005, "P", "AOSample");

	HW2_Test1(S, EOSample, file);
	HW2_Test2(S, AOSample, file);

	file << endl;
}

