//	Source file for Option class

#include "Option.hpp"

using namespace std;


//	Parent Option class!

//	Default constructor.
Option::Option() {}

//	Constructor.
Option::Option(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type, string new_name) :
	K(new_K), T(new_T), r(new_r), sig(new_sig), q(new_q), name(new_name) 
{
	if (new_type == "C" || new_type == "c" || new_type == "CALL" || new_type == "Call" || new_type == "call") type = "C";
	else if (new_type == "P" || new_type == "p" || new_type == "PUT" || new_type == "Put" || new_type == "put") type = "P";
}

//	Copy constructor.
Option::Option(const Option& o) :
	K(o.K), T(o.T), r(o.r), sig(o.sig), q(o.q), name(o.name)
{
	if (o.type == "C" || o.type == "c" || o.type == "CALL" || o.type == "Call" || o.type == "call") type = "C";
	else if (o.type == "P" || o.type == "p" || o.type == "PUT" || o.type == "Put" || o.type == "put") type = "P";
}

//	Destructor.
Option::~Option() {}

//	Assignment Operator.
Option& Option::operator = (const Option& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Get K.
double Option::Strike()
{
	return K;
}

//	Set K.
void Option::Strike(double new_K)
{
	K = new_K;
}

//	Get T.
double Option::Maturity()
{
	return T;
}

//	Set T.
void Option::Maturity(double new_T)
{
	T = new_T;
}

//	Get r.
double Option::InterestRate()
{
	return r;
}

//	Set r.
void Option::InterestRate(double new_r)
{
	r = new_r;
}

//	Get sig.
double Option::Volatility()
{
	return sig;
}

//	Set sig.
void Option::Volatility(double new_sig)
{
	sig = new_sig;
}

//	Get q.
double Option::Dividend()
{
	return q;
}

//	Set q.
void Option::Dividend(double new_q)
{
	q = new_q;
}

//	Get name.
string Option::Name()
{
	return name;
}

//	Set name.
void Option::Name(string new_name)
{
	name = new_name;
}

//	Get type.
string Option::Type()
{
	return type;
}

//	Set type.
void Option::Type(string new_type)
{
	if (new_type == "C" || new_type == "c" || new_type == "CALL" || new_type == "Call" || new_type == "call") type = "C";
	else if (new_type == "P" || new_type == "p" || new_type == "PUT" || new_type == "Put" || new_type == "put") type = "P";
}

//	Switch type.
void Option::Toggle()
{
	if (type == "C") type = "P";
	else type = "C";
}



//	European Options!

//	Default constructor.
EuropeanOption::EuropeanOption() {}

//	Constructor.
EuropeanOption::EuropeanOption(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type, string new_name) :
	Option(new_K, new_T, new_r, new_sig, new_q, new_type, new_name) {}

//	Copy constructor.
EuropeanOption::EuropeanOption(const EuropeanOption& o) :
	Option(o.K, o.T, o.r, o.sig, o.q, o.type, o.name) {}

//	Destructor.
EuropeanOption::~EuropeanOption() {}

//	Assignment Operator.
EuropeanOption& EuropeanOption::operator = (const EuropeanOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Calculate price with Black-Scholes-Merton.
double EuropeanOption::BSPrice(double S)
{
	if (type == "C") return BSCallPrice(S, K, T, r, sig, q);
	else return BSPutPrice(S, K, T, r, sig, q);
}

//	Calculate Delta with Black-Scholes-Merton.
double EuropeanOption::BSDelta(double S)
{
	if (type == "C") return BSCallDelta(S, K, T, r, sig, q);
	else return BSPutDelta(S, K, T, r, sig, q);
}

//	Calculate Gamma with Black-Scholes-Merton.
double EuropeanOption::BSGamma(double S)
{
	if (type == "C") return BSCallGamma(S, K, T, r, sig, q);
	else return BSPutGamma(S, K, T, r, sig, q);
}

//	Calculate Theta with Black-Scholes-Merton.
double EuropeanOption::BSTheta(double S)
{
	if (type == "C") return BSCallTheta(S, K, T, r, sig, q);
	else return BSPutTheta(S, K, T, r, sig, q);
}

//	European Option Black-Scholes Pricer.
vector<double> EuropeanOption::BSPricer(double S)
{
	return EuropeanBSPricer(type, S, K, T, r, sig, q);
}

//	Calculate Implied Volatility with Black-Scholes-Merton.
double EuropeanOption::BSImpVol(double V, double S, double acc)
{
	if (type == "C") return BSCallImpVol(V, S, K, T, r, sig, q, acc);
	else return BSPutImpVol(V, S, K, T, r, sig, q, acc);

}


//	Binomial Tree Pricer for European Options.
vector<double> EuropeanOption::BinomialTreePricer(double S, int N)
{
	return EuropeanBinomialTreePricer(type, S, K, T, r, sig, q, N);
}

//	Average Binomial Tree Pricer for European Options.
vector<double> EuropeanOption::ABTPricer(double S, int N)
{
	return EuropeanABTPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanOption::BBSPricer(double S, int N)
{
	return EuropeanBBSPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanOption::BBSRPricer(double S, int N)
{
	return EuropeanBBSRPricer(type, S, K, T, r, sig, q, N);
}


//	Trinomial Tree Pricer for European Options.
vector<double> EuropeanOption::TrinomialTreePricer(double S, int N)
{
	return EuropeanTrinomialTreePricer(type, S, K, T, r, sig, q, N);
}

//	Average Trinomial Tree Pricer for European Options.
vector<double> EuropeanOption::ATTPricer(double S, int N)
{
	return EuropeanATTPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes Pricer for European Options.
vector<double> EuropeanOption::TBSPricer(double S, int N)
{
	return EuropeanTBSPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for European Options.
vector<double> EuropeanOption::TBSRPricer(double S, int N)
{
	return EuropeanTBSRPricer(type, S, K, T, r, sig, q, N);
}


//	Monte-Carlo Pricer for European Options.
vector<double> EuropeanOption::MonteCarloPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanMonteCarloPricer(type, S, K, T, r, sig, q, N, Generator, seed);
}

//	Monte-Carlo Pricer for European Options.	//	Control Variate Technique.
double EuropeanOption::MCCVPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanMCCVPricer(type, S, K, T, r, sig, q, N, Generator, seed);
}

//	Monte-Carlo Pricer for European Options.	//	Antithetic Variablews.
double EuropeanOption::MCAVPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanMCAVPricer(type, S, K, T, r, sig, q, N, Generator, seed);
}

//	Monte-Carlo Pricer for European Options.	//	Moment Matching.
double EuropeanOption::MCMMPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanMCMMPricer(type, S, K, T, r, sig, q, N, Generator, seed);
}

//	Monte-Carlo Pricer for European Options.	//	Simultaneous Moment Matching and Control Variate Technique.
double EuropeanOption::MCCVMMPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanMCCVMMPricer(type, S, K, T, r, sig, q, N, Generator, seed);
}



//	American Options!

//	Default constructor.
AmericanOption::AmericanOption() {}

//	Constructor.
AmericanOption::AmericanOption(double new_K, double new_T, double new_r, double new_sig, double new_q, string new_type, string new_name) :
	Option(new_K, new_T, new_r, new_sig, new_q, new_type, new_name) {}

//	Copy constructor.
AmericanOption::AmericanOption(const AmericanOption& o) :
	Option(o.K, o.T, o.r, o.sig, o.q, o.type, o.name) {}

//	Destructor.
AmericanOption::~AmericanOption() {}

//	Assignment Operator.
AmericanOption& AmericanOption::operator = (const AmericanOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Binomial Tree Pricer for American Options.
vector<double> AmericanOption::BinomialTreePricer(double S, int N)
{
	return AmericanBinomialTreePricer(type, S, K, T, r, sig, q, N);
}

//	Average Binomial Tree Pricer for American Options.
vector<double> AmericanOption::ABTPricer(double S, int N)
{
	return AmericanABTPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanOption::BBSPricer(double S, int N)
{
	return AmericanBBSPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanOption::BBSRPricer(double S, int N)
{
	return AmericanBBSRPricer(type, S, K, T, r, sig, q, N);
}


//	Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::BTVRPricer(double S, int N)
{
	return AmericanBTVRPricer(type, S, K, T, r, sig, q, N);
}

//	Average Binomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::ABTVRPricer(double S, int N)
{
	return AmericanABTVRPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::BBSVRPricer(double S, int N)
{
	return AmericanBBSVRPricer(type, S, K, T, r, sig, q, N);
}

//	Binomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::BBSRVRPricer(double S, int N)
{
	return AmericanBBSRVRPricer(type, S, K, T, r, sig, q, N);
}


//	Trinomial Tree Pricer for American Options.
vector<double> AmericanOption::TrinomialTreePricer(double S, int N)
{
	return AmericanTrinomialTreePricer(type, S, K, T, r, sig, q, N);
}

//	Average Trinomial Tree Pricer for American Options.
vector<double> AmericanOption::ATTPricer(double S, int N)
{
	return AmericanATTPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes Pricer for American Options.
vector<double> AmericanOption::TBSPricer(double S, int N)
{
	return AmericanTBSPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.
vector<double> AmericanOption::TBSRPricer(double S, int N)
{
	return AmericanTBSRPricer(type, S, K, T, r, sig, q, N);
}


//	Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::TTVRPricer(double S, int N)
{
	return AmericanTTVRPricer(type, S, K, T, r, sig, q, N);
}

//	Average Trinomial Tree Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::ATTVRPricer(double S, int N)
{
	return AmericanATTVRPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::TBSVRPricer(double S, int N)
{
	return AmericanTBSVRPricer(type, S, K, T, r, sig, q, N);
}

//	Trinomial Tree Black-Scholes with Richardson Extrapolation Pricer for American Options.	//	With variance reduction.
vector<double> AmericanOption::TBSRVRPricer(double S, int N)
{
	return AmericanTBSRVRPricer(type, S, K, T, r, sig, q, N);
}



//	Barrier Options!

//	Default constructor.
BarrierOption::BarrierOption() {}

//	Constructor.
BarrierOption::BarrierOption(double new_K, double new_T, double new_r, double new_sig, double new_q, double new_B, string new_type, string new_Btype, string new_name) :
	Option(new_K, new_T, new_r, new_sig, new_q, new_type, new_name), B(new_B)
{
	if (new_Btype == "uni" || new_Btype == "UnI" || new_Btype == "UNI" || new_Btype == "uai" || new_Btype == "UaI" || new_Btype == "UAI" || new_Btype == "up and in" || new_Btype == "up-and-in" || new_Btype == "Up and In" || new_Btype == "Up-and-In" || new_Btype == "Up And In" || new_Btype == "Up-And-In" || new_Btype == "UP AND IN" || new_Btype == "UP-AND-IN")
		Btype = "UnI"; 
	else if (new_Btype == "dni" || new_Btype == "DnI" || new_Btype == "DNI" || new_Btype == "dai" || new_Btype == "DaI" || new_Btype == "DAI" || new_Btype == "down and in" || new_Btype == "down-and-in" || new_Btype == "Down and In" || new_Btype == "Down-and-In" || new_Btype == "Down And In" || new_Btype == "Down-And-In" || new_Btype == "DOWN AND IN" || new_Btype == "DOWN-AND-IN")
		Btype = "DnI";
	else if (new_Btype == "uno" || new_Btype == "UnO" || new_Btype == "UNO" || new_Btype == "uao" || new_Btype == "UaO" || new_Btype == "UAO" || new_Btype == "up and out" || new_Btype == "up-and-out" || new_Btype == "Up and Out" || new_Btype == "Up-and-Out" || new_Btype == "Up And Out" || new_Btype == "Up-And-Out" || new_Btype == "UP AND OUT" || new_Btype == "UP-AND-OUT")
		Btype = "UnO";
	else if (new_Btype == "dno" || new_Btype == "DnO" || new_Btype == "DNO" || new_Btype == "dao" || new_Btype == "DaO" || new_Btype == "DAO" || new_Btype == "down and out" || new_Btype == "down-and-out" || new_Btype == "Down and Out" || new_Btype == "Down-and-Out" || new_Btype == "Down And Out" || new_Btype == "Down-And-Out" || new_Btype == "DOWN AND OUT" || new_Btype == "DOWN-AND-OUT")
		Btype = "DnO";
}

//	Copy constructor.
BarrierOption::BarrierOption(const BarrierOption& o) :
	Option(o.K, o.T, o.r, o.sig, o.q, o.type, o.name), B(o.B), Btype(o.Btype) {}

//	Destructor.
BarrierOption::~BarrierOption() {}

//	Assignment Operator.
BarrierOption& BarrierOption::operator = (const BarrierOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		B = o.B;
		type = o.type;
		Btype = o.Btype;
		name = o.name;
	}
	return *this;
}


//	Get B.
double BarrierOption::Barrier()
{
	return B;
}

//	Set B.
void BarrierOption::Barrier(double new_B)
{
	B = new_B;
}

//	Get Btype.
string BarrierOption::BarrierType()
{
	return Btype;
}

//	Set Btype.
void BarrierOption::BarrierType(string new_Btype)
{
	Btype = new_Btype;
}


//	Calculate price with Black-Scholes-Merton.
double BarrierOption::BSPrice(double S)
{
	if (type == "C")
	{
		if (Btype == "DnO") return BSDnOCallPrice(S, K, T, r, sig, q, B);
	}
	return 0;	//			other types need to be added!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}


//	Monte-Carlo Pricer for Barrier Options.
double BarrierOption::MonteCarloPricer(double S, long int n, long int m, string Generator, long int seed)
{
	return BarrierMonteCarloPricer(type, S, K, T, r, sig, q, B, Btype, n, m, Generator, seed);
}

//	Monte-Carlo Pricer for Barrier Options.	//	With optimal path number and time interval number.
double BarrierOption::MonteCarloPricer(double S, long int N, string Generator, long int seed)
{
	return BarrierMonteCarloPricer(type, S, K, T, r, sig, q, B, Btype, N, Generator, seed);

}



//	Basket Options!

//	Default constructor.
BasketOption::BasketOption() {}

//	Constructor.
BasketOption::BasketOption(double new_K, double new_T, double new_r, double new_sig1, double new_sig2, double new_q1, double new_q2, double new_rho, string new_type,  string new_name) :
	Option(new_K, new_T, new_r, new_sig1, new_q1, new_type, new_name), sig2(new_sig2), q2(new_q2), rho(new_rho) {}

//	Copy constructor.
BasketOption::BasketOption(const BasketOption& o) :
	Option(o.K, o.T, o.r, o.sig, o.q, o.type, o.name), sig2(o.sig2), q2(o.q2), rho(o.rho) {}

//	Destructor.
BasketOption::~BasketOption() {}

//	Assignment Operator.
BasketOption& BasketOption::operator = (const BasketOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		sig2 = o.sig2;
		q2 = o.q2;
		rho = o.rho;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Get sig.
double BasketOption::Volatility(bool n)
{
	if (n == false) return sig;
	else return sig2;
}

//	Set sig.
void BasketOption::Volatility(double new_sig, bool n)
{
	if (n == false) sig = new_sig;
	else sig2 = new_sig;
}

//	Get q.
double BasketOption::Dividend(bool n)
{
	if (n == false) return q;
	else return q2;
}

//	Set q.
void BasketOption::Dividend(double new_q, bool n)
{
	if (n == false) q = new_q;
	else q2 = new_q;
}

//	Get rho.
double BasketOption::Rho()
{
	return rho;
}

//	Set rho.
void BasketOption::Rho(double new_rho)
{
	rho = new_rho;
}


//	Monte-Carlo Pricer for Basket Options.
double BasketOption::MonteCarloPricer(double S1, double S2, long int N, string Generator, long int seed)
{
	return BasketMonteCarloPricer(type, S1, S2, K, T, r, sig, sig2, rho, q, q2, N, Generator, seed);
}



//	Lookback Basket Options!

//	Default constructor.
LookbackBasketOption::LookbackBasketOption() {}

//	Constructor.
LookbackBasketOption::LookbackBasketOption(double new_K, double new_T, double new_r, double new_sig1, double new_sig2, double new_q1, double new_q2, double new_rho, string new_type, string new_name) :
	BasketOption(new_K, new_T, new_r, new_sig1, new_sig2, new_q1, new_q2, new_rho, new_type, new_name) {}

//	Copy constructor.
LookbackBasketOption::LookbackBasketOption(const LookbackBasketOption& o) :
	BasketOption(o.K, o.T, o.r, o.sig, o.sig2, o.q, o.q2, o.rho, o.type, o.name) {}

//	Destructor.
LookbackBasketOption::~LookbackBasketOption() {}

//	Assignment Operator.
LookbackBasketOption& LookbackBasketOption::operator = (const LookbackBasketOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		q = o.q;
		sig2 = o.sig2;
		q2 = o.q2;
		rho = o.rho;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Monte-Carlo Pricer for Lookback Basket Options.
double LookbackBasketOption::MonteCarloPricer(double S1, double S2, long int n, long int m, string Generator, long int seed)
{
	return LookbackBasketMonteCarloPricer(type, S1, S2, K, T, r, sig, sig2, rho, q, q2, n, m, Generator, seed);
}

//	Monte-Carlo Pricer for Lookback Basket Options.	//	With optimal path number and time interval number.
double LookbackBasketOption::MonteCarloPricer(double S1, double S2, long int N, string Generator, long int seed)
{
	return LookbackBasketMonteCarloPricer(type, S1, S2, K, T, r, sig, sig2, rho, q, q2, N, Generator, seed);
}



//	European Options with stock price follows Heston Model!

//	Default constructor.
EuropeanHestonOption::EuropeanHestonOption() {}

//	Constructor.
EuropeanHestonOption::EuropeanHestonOption(double new_K, double new_T, double new_r, double new_sig, double new_lambda, double new_sigbar, double new_eta, double new_rho, double new_q, string new_type, string new_name) :
	EuropeanOption(new_K, new_T, new_r, new_sig, new_q, new_type, new_name), lambda(new_lambda), sigbar(new_sigbar), eta(new_eta), rho(new_rho) {}

//	Copy constructor.
EuropeanHestonOption::EuropeanHestonOption(const EuropeanHestonOption& o) :
	EuropeanOption(o.K, o.T, o.r, o.sig, o.q, o.type, o.name), lambda(o.lambda), sigbar(o.sigbar), eta(o.eta), rho(o.rho) {}

//	Destructor.
EuropeanHestonOption::~EuropeanHestonOption() {}

//	Assignment Operator.
EuropeanHestonOption& EuropeanHestonOption::operator = (const EuropeanHestonOption& o)
{
	if (this != &o)
	{
		K = o.K;
		T = o.T;
		r = o.r;
		sig = o.sig;
		lambda = o.lambda;
		sigbar = o.sigbar;
		eta = o.eta;
		rho = o.rho;
		q = o.q;
		type = o.type;
		name = o.name;
	}
	return *this;
}


//	Get lambda.
double EuropeanHestonOption::Lambda()
{
	return lambda;
}

//	Set lambda.
void EuropeanHestonOption::Lambda(double new_lambda)
{
	lambda = new_lambda;
}

//	Get sigbar.
double EuropeanHestonOption::MeanVolatility()
{
	return sigbar;
}

//	Set sigbar.
void EuropeanHestonOption::MeanVolatility(double new_sigbar)
{
	sigbar = new_sigbar;
}

//	Get eta.
double EuropeanHestonOption::Eta()
{
	return eta;
}

//	Set eta.
void EuropeanHestonOption::Eta(double new_eta)
{
	eta = new_eta;
}

//	Get rho.
double EuropeanHestonOption::Rho()
{
	return rho;
}

//	Set rho.
void EuropeanHestonOption::Rho(double new_rho)
{
	rho = new_rho;
}


//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.
double EuropeanHestonOption::MonteCarloPricer(double S, long int n, long int m, string Generator, long int seed)
{
	return EuropeanHestonMonteCarloPricer(type, S, K, T, r, sig, lambda, sigbar, eta, rho, q, n, m, Generator, seed);
}

//	Monte-Carlo Pricer for European Options.	//	Stock price follows Heston Model.	//	With optimal path number and time interval number.
double EuropeanHestonOption::MonteCarloPricer(double S, long int N, string Generator, long int seed)
{
	return EuropeanHestonMonteCarloPricer(type, S, K, T, r, sig, lambda, sigbar, eta, rho, q, N, Generator, seed);
}


