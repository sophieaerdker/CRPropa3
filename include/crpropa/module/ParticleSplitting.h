#ifndef CRPROPA_PARTICLESPLITTING_H
#define CRPROPA_PARTICLESPLITTING_H

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "kiss/logger.h"


namespace crpropa {

/**
@class ParticleSplitting
@brief Candidates are split into n copies when they cross specified energy bins. Weights are set accordingly.
*/

class ParticleSplittingModule: public Module {
private:
	double n_split;
	std::vector<double> Ebins;

public:
	
	ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param n_bins		Number of energy bins 
	 */
	ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins, bool log);
	/** Constructor
	 @param n_split 	Number of copies candidates are split 
	 @param Emin 		Minimal energy for splitting
	 @param Emax		Maximal energy for splitting
	 @param n_bins		Number of energy bins 
	 @param log 		Energy bins in log
	 */
	void process(Candidate *c) const;

	void setEnergyBins(double Emin, double Emax, double n_bins, bool log);

	void setNsplit(int n);

};
/** @}*/



}; // end namesspace crpropa
#endif // CRPROPA_PARTICLESPLITTING_H
