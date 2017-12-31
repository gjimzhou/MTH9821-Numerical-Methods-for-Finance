//	Main program for testing

#include <iostream>
#include <iomanip> 
#include <limits>

#include "BinomialTree.hpp"


int main()
{
	int precision = numeric_limits<double>::max_digits10;
	cout << setprecision(precision) << endl;

	double S = 50;
	double K = 55.55;
	double T1 = 3.0 / 12.0;
	double T2 = 7.0 / 12.0;
	double r = 0.02;
	double sig = 0.30;

	int N1 = 6;
	int N2 = 7;
	double tol = 0.0001;

	vector<Dividend> q0;

	vector<Dividend> q1;
	q1.push_back(Dividend(0.99, 2.0 / 12.0, "p"));
	q1.push_back(Dividend(0.99, 4.0 / 12.0, "p"));
	q1.push_back(Dividend(0.99, 6.0 / 12.0, "p"));

	vector<Dividend> q2;
	q2.push_back(Dividend(0.50, 2.0 / 12.0, "f"));
	q2.push_back(Dividend(0.50, 4.0 / 12.0, "f"));
	q2.push_back(Dividend(0.50, 6.0 / 12.0, "f"));

	vector<Dividend> q3;
	q3.push_back(Dividend(0.50, 2.0 / 12.0, "f"));
	q3.push_back(Dividend(0.99, 4.0 / 12.0, "p"));
	q3.push_back(Dividend(0.72, 6.0 / 12.0, "f"));

	vector<double> res;
	double value_old;
	double value_new;
	double diff;
	int N;

	//	3M EP F.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N1;
	
	while (diff > tol)
	{
		N *= 2;
		vector<double> V = EuropeanDividendBinomialTreePricer("P", S, K, T1, r, sig, q1, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	3M AP F.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N1;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = AmericanDividendBinomialTreePricer("P", S, K, T1, r, sig, q1, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	7M EP F.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = EuropeanDividendBinomialTreePricer("P", S, K, T2, r, sig, q1, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	7M AP F.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = AmericanDividendBinomialTreePricer("P", S, K, T2, r, sig, q1, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;


	//	3M EP P.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N1;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = EuropeanDividendBinomialTreePricer("P", S, K, T1, r, sig, q2, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	3M AP P.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N1;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = AmericanDividendBinomialTreePricer("P", S, K, T1, r, sig, q2, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	7M EP P.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = EuropeanDividendBinomialTreePricer("P", S, K, T2, r, sig, q2, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	7M AP P.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = AmericanDividendBinomialTreePricer("P", S, K, T2, r, sig, q2, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;


	//	7M EP M.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = EuropeanDividendBinomialTreePricer("P", S, K, T2, r, sig, q3, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;

	//	7M AP M.
	value_old = 0;
	value_new = 1;
	diff = 1;
	N = N2;

	while (diff > tol)
	{
		N *= 2;
		vector<double> V = AmericanDividendBinomialTreePricer("P", S, K, T2, r, sig, q3, N);
		res = V;
		value_old = value_new;
		value_new = V[0];
		diff = abs(value_new - value_old);
	}
	cout << N << " ";
	for (auto& v : res)
	{
		cout << v << " ";
	}
	cout << endl;



	//double SNon = DivToNon(q1, S, T1, r);
	//double SDiv = NonToDiv(q1, SNon, T1, r);
	//cout << S << " " << SNon << " " << SDiv << endl;

	system("pause");
	return 0;
}