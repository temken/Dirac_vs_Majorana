#ifndef __Targets_hpp_
#define __Targets_hpp_

#include <vector>
#include <string>

struct Atomic_Electron 
{
	unsigned int n,l;

	// Electronic response tables
	double k_min,k_max,q_min,q_max;
	double logk_min, logk_max,logq_min,logq_max, dlogk, dlogq;
	unsigned int Nk, Nq;
	std::vector<double> k_Grid={};
	std::vector<double> q_Grid={};
	std::string Name;
	double Binding_Energy;
	double Nucleus_Mass;
	int ne_Secondary;

	Atomic_Electron(std::string element,double A, int N, int L, double Ebinding, double kMin,double kMax, double qMin, double qMax,int neSecondary = 0);

	std::vector< std::vector< std::vector<double> > > Response_Tables;
	double Electronic_Response(unsigned int response,double q, double E) const;
};

struct Atom
{
	std::string Name;

	int Z;
	double A;
	
	double mass;

	double W;

	std::vector<Atomic_Electron> electrons;

	//Constructor
	Atom()
	:Name("default"), Z(0), A(0), mass(0), W(0), electrons({})
	{};
	Atom(std::string element_name, int z, double a, std::vector<Atomic_Electron> shells ={});

	double Lowest_Binding_Energy() const;

	Atomic_Electron Electron(int n, int l);

};

extern Atom Import_Electronic_Responses(std::string element);

#endif