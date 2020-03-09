#ifndef __Direct_Detection_hpp_
#define __Direct_Detection_hpp_

#include <iostream>
#include <string>

#include "Numerics_Functions.hpp"
#include "DM_Particle.hpp"	


//1. Detector base class
class Detector
{
	protected:
		
		double flat_efficiency;
		long unsigned int observed_signals;
		double certainty_level;

		virtual double Maximum_Energy_Deposit(const DM_Particle& DM) const { return 0.0; };
		virtual double Minimum_DM_Mass(DM_Particle& DM) const {return 0.0;};
		Interpolation Spectrum_Base(const DM_Particle& DM,double Emin,double Emax,unsigned int points = 200);
		void Print_Summary_Base() const;
		
	public:
		std::string name;
		double exposure; 
		std::string target_particles;
		Detector() :  flat_efficiency(1.0), observed_signals(0), certainty_level(0.9),name("Generic"), exposure(0.0), target_particles("default") {};
		Detector(double expo,std::string target_type) : flat_efficiency(1.0), observed_signals(0), certainty_level(0.9), name("Generic"), exposure(expo) ,target_particles(target_type) {};

		void Set_Name(std::string n);
		void Set_Flat_Efficiency(double eff);
		void Set_Observed_Signals(unsigned long int n);
		void Set_Certainty_Level(double cl);

		virtual std::vector<std::vector<double>> Limit_Curve(DM_Particle& DM,double mMin,double mMax, int points);

		virtual double dRdE(double E, const DM_Particle& DM,double vDM = 1e-3) { return 0.0;};
		virtual Interpolation Spectrum(const DM_Particle& DM,unsigned int points = 200) {return Interpolation();};
		
		virtual double N_Signals(const DM_Particle& DM) { return 0.0; };
		virtual double Likelihood(const DM_Particle& DM) { return 0.0;};
		virtual double Upper_Bound(DM_Particle& DM) {return 0.0;};

		virtual double Minimum_Speed(const DM_Particle& DM) const {return 0.0;};
		
		virtual void Print_Summary() const {Print_Summary_Base();};
};

#endif