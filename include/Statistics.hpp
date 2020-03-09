#ifndef __Statistics_hpp_
#define __Statistics_hpp_

#include <vector>
#include <random>

//1. Distributions
extern double PDF_Uniform(double x, double x_min, double x_max);
extern double CDF_Uniform(double x, double x_min, double x_max);

extern double PDF_Gauss(double x, double mu, double sigma);
extern double CDF_Gauss(double x, double mu, double sigma);
extern double Quantile_Gauss(double p, double mu, double sigma);

extern double PMF_Binomial(int mu, double p, int trials);
extern double CDF_Binomial();

extern double PMF_Poisson(double expected_events, unsigned int events);
extern double CDF_Poisson(double expectation_value,unsigned int observed_events);
extern double Inv_CDF_Poisson(unsigned int observed_events, double cdf); //Solves the CDF = cdf for mu

extern double PDF_Chi_Square(double x, double dof);
extern double CDF_Chi_Square(double x, double dof);
extern double PDF_Chi_Bar_Square(double x, std::vector<double> weights);
extern double CDF_Chi_Bar_Square(double x, std::vector<double> weights);

//2. Sampling random numbers
extern double Sample_Uniform(std::mt19937& PRNG, double x_min = 0.0, double x_max = 1.0);

extern int Sample_Poisson(std::mt19937& PRNG, double expectation_value);

//3. Likelihoods
extern double Likelihood_Poisson_Binned(std::vector<double> expectation_values,std::vector<int> observed_events);
extern double Log_Likelihood_Poisson_Binned(std::vector<double> expectation_values,std::vector<int> observed_events);

//4. Basic data analysis
extern double Arithmetic_Mean(const std::vector<double>& data);
extern double Median(std::vector<double>& data);
extern double Variance(const std::vector<double>& data);
extern double Standard_Deviation(const std::vector<double>& data);


#endif