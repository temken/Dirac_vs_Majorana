#ifndef __Direct_Detection_Electron_hpp_
#define __Direct_Detection_Electron_hpp_

#include <string>

#include "DM_Particle.hpp"
#include "Targets.hpp"
#include "Direct_Detection.hpp"

extern double vMinimal_e(double q,double Ee,double mDM);
extern double dR_dlnE_Ion(double E,const DM_Particle& DM, const Atomic_Electron& shell);

extern double dR_dlnE_Ion(double E,const DM_Particle& DM, const Atom& atom);
extern double dR_dlnE_Ion(unsigned int response,double E,const DM_Particle& DM, const Atom& atom);

extern double R_Total(const DM_Particle& DM, const Atom& atom);


class Detector_Ionization : public Detector
{
	private:
		Atom target;

		unsigned long int observed_signals;

		bool binned_data;
		std::vector<double> bin_efficiencies;
		std::vector<unsigned long int> binned_events;

		virtual double Maximum_Energy_Deposit(const DM_Particle& DM) const override;
		virtual double Minimum_DM_Mass(DM_Particle& DM) const override;

		bool Using_Electron_Bins;
		
		double PDFne(unsigned int ne,double Ee,const Atomic_Electron& shell) const;

		bool Using_PE_Bins;
		double mu_PE, sigma_PE;
		std::vector<double> Trigger_Efficiency_PE;
		std::vector<double> Acceptance_Efficiency_PE;
		std::vector<int> bins;

	public:
		unsigned int ne_threshold;
		
		Detector_Ionization(std::string atom,double expo);

		Detector_Ionization(double expo);

		void Import_Atomic_Response_Functions(std::string atom);
		
		virtual double dRdE(double E, const DM_Particle& DM,double vDM = 1e-3) override;
		virtual Interpolation Spectrum(const DM_Particle& DM,unsigned int points = 200) override;

		double dRdne(unsigned int ne,const DM_Particle& DM,const Atomic_Electron& shell);
		// double dRdne_2(unsigned int ne,const DM_Particle& DM,const Atomic_Electron& shell);
		double dRdne(unsigned int ne,const DM_Particle& DM);
		
		void Set_Electron_Threshold(unsigned int n);
		void Set_Binned_Events(const std::vector<unsigned long int>& events,std::vector<double> eff = {});


		void Set_PE_Distribution(double mu,double sigma, std::vector<int> binsizes);
		void Import_Trigger_Efficiency_PE(std::string filename);
		void Import_Acceptance_Efficiency_PE(std::string filename);
		double dRdS2(unsigned int S2,const DM_Particle& DM);

		
		virtual double N_Signals(const DM_Particle& DM) override;
		virtual double Likelihood(const DM_Particle& DM) override;
		virtual double Upper_Bound(DM_Particle& DM) override;

		virtual std::vector<std::vector<double>> Limit_Curve(DM_Particle& DM,double mMin,double mMax, int points) override;


		virtual double Minimum_Speed(const DM_Particle& DM) const override;
		
		virtual void Print_Summary() const override;

};

extern void test_function(int argument);
#endif