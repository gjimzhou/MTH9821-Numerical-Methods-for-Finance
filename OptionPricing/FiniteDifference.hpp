//	Header file for Finite Difference Pricers

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "NumericalMethod.hpp"
#include "OptionFunc.hpp"

using namespace std;


#ifndef FiniteDifference_HPP
#define FiniteDifference_HPP

//	Finite Difference Pricer for European Options.
vector<double> EuropeanFiniteDifferencePricer(string type, double S, double K, double T, double r, double sig, double q, long int M, double alpha);



#endif
