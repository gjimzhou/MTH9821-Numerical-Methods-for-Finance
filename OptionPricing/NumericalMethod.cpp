//	Source file for Numerical methods functions

#include "NumericalMethod.hpp"

using namespace std;

const double PI = 3.141592653589793238463;


//	CDF for standard normal distribution.
double pnorm(double x)
{
	return 0.5 * (1 + erf(x / sqrt(2)));
}

//	PDF for standard normal distribution.
double dnorm(double x)
{
	return exp(-x * x / 2) / sqrt(2 * PI);
}

//	QF for standard normal distribution.
double qnorm(double u)
{
	double a0 = 2.50662823884;
	double a1 = -18.61500062529;
	double a2 = 41.39119773534;
	double a3 = -25.44106049637;

	double b0 = -8.47351093090;
	double b1 = 23.08336743743;
	double b2 = -21.06224101826;
	double b3 = 3.13082909833;

	double c0 = 0.3374754822726147;
	double c1 = 0.9761690190917186;
	double c2 = 0.1607979714918209;
	double c3 = 0.0276438810333863;
	double c4 = 0.0038405729373609;
	double c5 = 0.0003951896511919;
	double c6 = 0.0000321767881768;
	double c7 = 0.0000002888167364;
	double c8 = 0.0000003960315187;

	double x;
	double y = u - 0.5;

	if (abs(y) < 0.42)
	{
		double r = y * y;
		x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1);
	}
	else
	{
		double r = u;
		if (y > 0) r = 1 - u;

		r = log(-log(r));
		x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));

		if (y < 0) x = -x;
	}

	return x;
}


//	Initialize random seed.
long int SeedInit(long int seed)
{
	if (seed == 0) seed = time(0);
	return seed;
}

//	Linear Congruential Generator for standard uniform distribution random variable.
vector<double> LinearCongruentialGenerator(long int N, long int seed)
{
	//	Value to begin iteration.
	seed = SeedInit(seed);

	long int m = 2147483647;	//	2^31 - 1.
	long int a = 39373;	//	Magic prime number.

	//	Method involve q and r proved at page 45-46 from Glasserman.
	long int q = m / a;
	long int r = m % a;

	vector<double> result;

	for (long int i = 0; i < N; i++)
	{
		long int k = seed / q;

		seed = a * (seed - k * q) - k * r;
		if (seed < 0) seed = seed + m;

		result.push_back(seed / (double)m);
	}

	return result;
}

// Inverse Transform Method for standard normal distribution random variable.
vector<double> InverseTransformMethod(long int N, long int seed)
{
	vector<double> runif = LinearCongruentialGenerator(N, seed);

	vector<double> result;

	for (long int i = 0; i < N; i++)
	{
		result.push_back(qnorm(runif[i]));
	}

	return result;
}

// Acceptance-Rejection Method for standard normal distribution random variable.
vector<double> AcceptanceRejectionMethod(long int N, long int seed)
{
	vector<double> result;

	while (result.size() < N)
	{
		long int M = result.size();

		//	Generate uniform to be used then reset the random seed.
		vector<double> runif = LinearCongruentialGenerator(3 * (N - M), seed);
		seed = runif[3 * (N - M) - 1] * 2147483647;

		for (long int i = 0; i < (N - M); i++)
		{
			double u1 = runif[3 * i];
			double u2 = runif[3 * i + 1];
			double u3 = runif[3 * i + 2];

			double x = -log(u1);

			if (u2 <= exp(-0.5 * (x - 1) * (x - 1)))
			{
				if (u3 <= 0.5) x = -x;

				result.push_back(x);
			}
		}
	}

	return result;
}

// Box-Muller Method for standard normal distribution random variable.
vector<double> BoxMullerMethod(long int N, long int seed)
{
	vector<double> result;

	while (result.size() < N)
	{
		long int M = result.size();

		//	Generate uniform to be used then reset the random seed.
		vector<double> runif = LinearCongruentialGenerator(N - M, seed);
		seed = runif[N - M - 1] * 2147483647;

		for (long int i = 0; i < (N - M); i += 2)
		{
			double u1 = 2 * runif[i] - 1;
			double u2 = 2 * runif[i + 1] - 1;

			double x = u1 * u1 + u2 * u2;

			if (x <= 1)
			{
				double y = sqrt(-2 * log(x) / x);

				result.push_back(u1 * y);
				result.push_back(u2 * y);
			}
		}
	}

	return result;
}


//	Generator for standard normal distribution random variables.
vector<double> rnorm(long int N, string Generator, long int seed)
{
	vector<double> z;

	if (Generator == "BMM")
	{
		z = BoxMullerMethod(N, seed);
	}
	else if (Generator == "ARM")
	{
		z = AcceptanceRejectionMethod(N, seed);
	}
	else
	{
		z = InverseTransformMethod(N, seed);
	}

	return z;
}


//	Forward Substitution.
vector<double> ForwardSubstitution(vector<vector<double>> L, vector<double> b)
{
	long int n = b.size();
	vector<double> x(n);

	x[0] = b[0] / L[0][0];
	for (long int j = 1; j < n; j++)
	{
		double sum = 0;
		for (long int k = 0; k < j; k++)
		{
			sum += L[j][k] * x[k];
		}
		x[j] = (b[j] - sum) / L[j][j];
	}

	return x;
}

//	Backward Substitution.
vector<double> BackwardSubstitution(vector<vector<double>> U, vector<double> b)
{
	long int n = b.size();
	vector<double> x(n);

	x[0] = b[n - 1] / U[n - 1][n - 1];
	for (long int j = n - 1; j > -1; j--)
	{
		double sum = 0;
		for (long int k = j + 1; k < n; k++)
		{
			sum += U[j][k] * x[k];
		}
		x[j] = (b[j] - sum) / U[j][j];
	}

	return x;
}

//	Forward Substitution for bidiagonal matrix.
vector<double> ForwardSubstitutionBidiag(vector<vector<double>> L, vector<double> b)
{
	long int n = b.size();
	vector<double> x(n);

	x[0] = b[0] / L[0][0];
	for (long int j = 1; j < n; j++)
	{
		x[j] = (b[j] - L[j][j - 1] * x[j - 1]) / L[j][j];
	}

	return x;
}

//	Backward Substitution for bidiagonal matrix.
vector<double> BackwardSubstitutionBidiag(vector<vector<double>> U, vector<double> b)
{
	long int n = b.size();
	vector<double> x(n);

	x[0] = b[n - 1] / U[n - 1][n - 1];
	for (long int j = n - 1; j > -1; j--)
	{
		x[j] = (b[j] - U[j][j + 1] * x[j + 1]) / U[j][j];
	}

	return x;
}