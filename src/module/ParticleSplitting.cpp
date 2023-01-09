#include "crpropa/module/ParticleSplitting.h"

namespace crpropa {

ParticleSplittingModule::ParticleSplittingModule() {
	// no particle splitting if EnergyBins and NSplit are not specified
	setNsplit(0);
}

ParticleSplittingModule::ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins) {
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, false);
}

ParticleSplittingModule::ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins, bool log) {
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, log);
}

ParticleSplittingModule::ParticleSplittingModule(int SpectralIndex, double Emin, int factor)  {
	// for use with Diffusive Shock Acceleration

	if (SpectralIndex < 0){
		throw std::runtime_error(
				"ParticleSplitting: spectralIndex < 0 !");
	}

	double Emax = Emin*pow(10, factor); 
	int n_bins = factor + 1;
	setEnergyBins(Emin, Emax, n_bins, true);
	// to compensate for loss of particles per energy bin:
	setNsplit((int) pow(10, SpectralIndex-1)); 

}


void ParticleSplittingModule::process(Candidate *c) const {

	double currE = c->current.getEnergy(); 
	double prevE = c->previous.getEnergy();

	if (currE < Ebins[0] || n_split == 0){
		// current energy is smaller than first bin -> no splitting
		// or, number of splits = 0
		return;
	}
	for (size_t i = 0; i < Ebins.size(); ++i){
		
		if( prevE < Ebins[i] ){
			// previous energy is in energy bin [i-1, i]
			if(currE < Ebins[i]){
				//assuming that dE greater than 0, prevE and E in same energy bin -> no splitting
				return;
			}

			// current energy is in energy bin [i,i+1] or higher -> particle splitting for each crossing
			for (size_t j = i; j < Ebins.size(); ++j ){

				// adapted from Acceleration Module:
				c->updateWeight(1. / n_split); // * 1/n_split

				for (int i = 1; i < n_split; i++) {
				
					ref_ptr<Candidate> new_candidate = c->clone(false);
					new_candidate->parent = c;
					uint64_t snr = Candidate::getNextSerialNumber();
					new_candidate->setSerialNumber(snr);
					c->addSecondary(new_candidate);
					Candidate::setNextSerialNumber(snr + 1);
					//std::cout<< "new serial number" << snr << std::endl;
				}


				if (j < Ebins.size()-1 && currE < Ebins[j+1]){
					// candidate is in energy bin [j, j+1] -> no further splitting
					return;
				}
			}

			return;

		}
	}
}
	

void ParticleSplittingModule::setEnergyBins(double Emin, double Emax, double n_bins, bool log) {

	Ebins.resize(0);

	if (Emin > Emax){
		throw std::runtime_error(
				"ParticleSplitting: Emin > Emax!");
	}

	double dE = (Emax-Emin)/n_bins;
	
	for (size_t i = 0; i < n_bins; ++i) {
		if (log == true) {
			Ebins.push_back(Emin * pow(Emax / Emin, i / (n_bins - 1.0)));
		} else {
			Ebins.push_back(Emin + i * dE);
		}
	}

	std::cout << Ebins[0] << Ebins[n_bins-1] << std::endl;
}

void ParticleSplittingModule::setNsplit(int n) {

	n_split = n;
}


} // end namespace crpropa

