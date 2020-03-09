#ifndef __Numerics_Functions_hpp_
#define __Numerics_Functions_hpp_

#include <functional>
#include <vector>

//1. Special functions
	//1.1 Simple functions
		extern int Sign(double arg);
		extern double Sign(double x, double y); //Returns x with the sign of y.
		extern double StepFunction(double x);
		extern double Round(double N,unsigned int digits=3);
	//1.2 Gamma Functions
		// extern std::vector<double> FactorialList;
		extern double Factorial(unsigned int n); 
		extern double Binomial_Coefficient(int n,int k);
		extern double GammaLn(double x);
		extern double Gamma(double x);
		extern double Upper_Incomplete_Gamma(double x, double s);
		extern double Lower_Incomplete_Gamma(double x, double s);
		extern double GammaQ(double x,double a); //Incomplete Gamma Q(x,a)
		extern double GammaP(double x,double a); //Incomplete Gamma P(x,a)
		extern double Inv_GammaP(double p,double a); //Solves P(x,a)=p for x.
		extern double Inv_GammaQ(double q,double a); //Solves Q(x,a)=q for x.

		extern double Inv_erf(double p);

//2. Numerical Integration via Adaptive Simpson Method
	extern double Find_Epsilon(std::function<double(double)> func, double a,double b,double precision);
	extern double Integrate(std::function<double(double)> func, double a,double b, double epsilon,int maxRecursionDepth=20);
	double auxAdaptiveSimpson(std::function<double(double)> func, double a,double b, double epsilon,double S,double fa,double fb,double fc,int bottom,bool &warning);

//3. Root finding
	extern double Find_Root(std::function<double(double)>& func,double xLeft, double xRight,double epsilon);

//4. One-dimensional interpolation of tabulated data using Steffen splines 
	class Interpolation 
	{
		private:
			std::vector<std::vector<double>> TabulatedData;
			unsigned int N_Data;
			std::vector<double> xDomain;
			//Steffen coefficients
			std::vector<double> a,b,c,d;
			//Pre-factor
			double preFactor;
			//compute Steffen coefficients
			void Compute_Steffen_Coefficients(std::vector<std::vector<double>>& data, std::vector<double> &a,std::vector<double> &b,std::vector<double> &c,std::vector<double> &d);
			//Locate j in array such that x[j]<x<x[j+1].
			unsigned int jLast;
			bool corr;// if successive calls are correlated, then the hunt method can be faster.
			unsigned int Bisection(double x,int jLeft, int jRight);
			unsigned int Hunt(double x);
			unsigned int Locate(double x);
		public:
			//Constructor from data or a data file
				Interpolation();
				Interpolation(const std::vector<std::vector<double>>& data,double dim1=1.0,double dim2=1.0);
				Interpolation(const std::string& filename,double dim1=1.0,double dim2=1.0);
			//Return values
				std::vector<std::vector<double>> Return_Data();
				std::vector<double> Return_Domain();
				double Return_Prefactor();
				std::vector<std::vector<double>> Return_Coefficients();
			//Set values
				void Set_Prefactor(double factor);
			//Interpolation
				double Interpolate(double x);
				double operator ()(double x)
			    {
			        return Interpolate(x);
			    };
			//Multiply by a constant
			    void Multiply(double factor);
			//Save function in a file
			    void Save_Function(std::string filename,unsigned int points);
	};
	
//5. Read in table into vector
	extern std::vector<double> Read_List(std::string filepath,double dimension=1.0);

//6. Save vector into file
	extern void Save_List(std::string filepath, std::vector<double> list, double dimension = 1.0);
	extern void Save_List(std::string filepath, const std::vector<std::vector<double>>& list, std::vector<double> dimensions ={});

//7. Create list with equi-distant numbers in log-space
	extern std::vector<double> Log_Space(double min, double max, unsigned int steps);


//8. Multidimensional minimization Amoeba
	struct Minimization
	{
		const double ftol;
		int nfunc; //The number of function evaluations.
		int mpts;
		int ndim;
		double fmin;	//Function value at the minimum.
		std::vector<double> y;	//Function values at the vertices of the simplex.
		std::vector<std::vector<double>> current_simplex;	// p
		Minimization(const double ftoll) : ftol(ftoll) {} //The constructor argument ftoll is the fractional convergence tolerance to be achieved in the function value (n.b.!).
		
		//Multidimensional minimization of the function or functor func(x), where x[0..ndim-1] is a vector in ndim dimensions, by the downhill simplex method of Nelder and Mead. The initial simplex is specified as in equation (10.5.1) by a point[0..ndim-1] and a constant displacement del along each coordinate direction. Returned is the location of the minimum.
		std::vector<double> minimize(std::vector<double> &starting_point, const double delta, std::function<double(std::vector<double>)> func);

		//Alternative interface that takes different displacements dels[0..ndim-1] in different di- rections for the initial simplex.
		std::vector<double> minimize(std::vector<double> &starting_point, std::vector<double> &deltas, std::function<double(std::vector<double>)> func);

		//Most general interface: initial simplex specified by the matrix pp[0..ndim][0..ndim-1]. Its ndim+1 rows are ndim-dimensional vectors that are the vertices of the starting simplex.
		std::vector<double> minimize(std::vector<std::vector<double>> &pp, std::function<double(std::vector<double>)> func);

		void get_psum(std::vector<std::vector<double>> &p, std::vector<double> &psum); //Utility function.

		double amotry(std::vector<std::vector<double>> &p, std::vector<double> &y, std::vector<double> &psum,const int ihi, const double fac, std::function<double(std::vector<double>)> func);
	};



#endif