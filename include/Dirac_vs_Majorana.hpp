#ifndef __Dirac_vs_Majorana_hpp_
#define __Dirac_vs_Majorana_hpp_

#include <string>
#include <random>

#include "Linear_Algebra.hpp"
#include "DM_Particle.hpp"
#include "Direct_Detection_Electron.hpp"

//1. Create a DM particle interacting with the photon via higher order operators (anapole, magnetic dipole, electric dipole)
extern DM_Particle Create_DM_Particle(double mDM, double C_1, double C_2, double C_3);

//2. Event spectrum in terms of n_e.
extern std::vector<double> Binned_Event_Spectrum(Detector_Ionization& detector, double mDM, double C_1, double C_2, double C_3);

//3. Tabulate the three fiducial spectra.
extern std::vector<std::vector<double>> Fiducial_Spectra(Detector_Ionization& detector, double mDM);

//4. Fix couplings to a given number of events and hierarchy.
extern std::vector<double> Determine_Couplings(Detector_Ionization& detector, double mDM, double number_of_events, std::vector<double> operator_hierarchy, std::vector<std::vector<double>> fiducial_spectra = {});

//5. Monte Carlo simulate an experiment generating number of events in a given number of bins between E_i, E_f
extern std::vector<std::vector<double>> Simulate_Experiment(std::mt19937& PRNG, Detector_Ionization& detector, double mDM, double C_1, double C_2, double C_3, std::vector<std::vector<double>> fiducial_spectra = {}, bool terminal_output = false);

extern std::vector<std::vector<double>> Import_MC_Data(std::string filename, bool terminal_output = false);

//6. Statistics functions
extern double Test_Statistic_t(std::vector<std::vector<double>>& data, Detector_Ionization& detector,double mDM, std::vector<std::vector<double>> fiducial_spectra = {}, bool terminal_output = false);

extern Matrix Fisher_Matrix_Asimov(Detector_Ionization& detector, double mDM,  double C_1, double C_2, double C_3, std::vector<std::vector<double>>& fiducial_spectra);

extern std::vector<double> Sample_t(unsigned int sample_size, std::mt19937& PRNG, Detector_Ionization& detector, double mDM, double C_1_true, double C_2_true, double C_3_true,  std::vector<std::vector<double>> fiducial_spectra = {}, bool terminal_output = false);

extern std::vector<double> Asymptotic_Distribution_Weights(std::vector<double> &t_sample, int N_weights);

extern double Compute_p_Value(unsigned int sample_size, std::mt19937& PRNG, Detector_Ionization& detector, double number_of_signals, double mDM, std::vector<double> operator_hierarchy,  std::vector<std::vector<double>> fiducial_spectra = {},  bool terminal_output = false);

extern double Significance(double p_value);

extern double Significance_Threshold(double Z, unsigned int sample_size, std::mt19937& PRNG, Detector_Ionization& detector, double mDM, std::vector<double> operator_hierarchy, double N_min = 1, double N_max = 1000, std::vector<std::vector<double>> fiducial_spectra = {},  bool terminal_output = false);

#endif