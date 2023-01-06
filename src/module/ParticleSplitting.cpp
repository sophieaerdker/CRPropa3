#include "crpropa/module/ParticleSplitting.h"

namespace crpropa {

ParticleSplittingModule::ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins) {
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, false);
}

ParticleSplittingModule::ParticleSplittingModule(int n_split, double Emin, double Emax, double n_bins, bool log) {
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, log);
}


void ParticleSplittingModule::process(Candidate *c) const {

	double currE = c->current.getEnergy(); 
	double prevE = c->previous.getEnergy();

	if (currE < Ebins[0]){
		// current energy is smaller than first bin -> no splitting
		return;
	}
	// in case previous energy is already in maxenergy bin, no splitting:
	for (size_t i = 1; i < Ebins.size()-1; ++i){
		
		if( prevE < Ebins[i] ){
			// previous energy is in energy bin [i-1, i]
			if(currE < Ebins[i]){
				//assuming that dE greater than 0, prevE and E in same energy bin -> no splitting
				//std::cout << "SAME BIN prevE: " << prevE << " , prev Ebin: " << Ebins[i] << " , currE: " << currE << " , curr Ebin: " << Ebins[i+1] << std::endl;
				return;
			}

			//std::cout << "BIN CROSSING prevE: " << prevE << " , prev Ebin: " << Ebins[i] << " , currE: " << currE << " , curr Ebin: " << Ebins[i+1] << std::endl;
			// current energy is in energy bin [i,i+1] or higher -> particle splitting for each crossing

			for (size_t j = i; j < Ebins.size()-1; ++j ){

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

				if (currE < Ebins[j+1]){
					// candidate is in energy bin [j, j+1] -> no further splitting
					return;
				}
			}	
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
}

void ParticleSplittingModule::setNsplit(int n) {

	n_split = n;
}


} // end namespace crpropa

