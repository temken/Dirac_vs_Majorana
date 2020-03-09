#ifndef __DM_Particle_hpp_
#define __DM_Particle_hpp_

#include "Physics_Functions.hpp"
#include "Targets.hpp"

//Base class for a DM particle with virtual functions for the cross sections
class DM_Particle
{
	protected:
		std::vector<double> couplings_h;
		std::vector<double> couplings_l;
	public:
		double mass;
		double spin;

		DM_Particle();
		DM_Particle(double mass, double spin=0.5);

		void Set_Mass(double m);
		void Set_Light_DM(bool ldm);

		void Set_Coupling(unsigned int index,double value_contact,double value_longrange = 0.0);
		void Reset_Couplings();

		double DM_Response(unsigned int response,double q, double dE,double vDM) const;
		double dSigmaV_dlnE(double E,const Atomic_Electron& shell) const;
		double dSigmaV_dlnE(unsigned int response, double E,const Atomic_Electron& shell) const;

		void Print_Summary() const;
};


#endif