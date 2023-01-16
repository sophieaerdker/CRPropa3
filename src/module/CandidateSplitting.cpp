#include "crpropa/module/CandidateSplitting.h"

namespace crpropa {

CandidateSplitting::CandidateSplitting() {
	// no particle splitting if EnergyBins and NSplit are not specified
	setNsplit(0);
}

CandidateSplitting::CandidateSplitting(int n_split, double Emin, double Emax, double n_bins) {
	// linear energy bins:
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, false);
}

CandidateSplitting::CandidateSplitting(int n_split, double Emin, double Emax, double n_bins, bool log) {
	setNsplit(n_split);
	setEnergyBins(Emin, Emax, n_bins, log);
}

CandidateSplitting::CandidateSplitting(int SpectralIndex, double Emin, int factor)  {
	// to use with Diffusive Shock Acceleration

	if (SpectralIndex <= 0){
		throw std::runtime_error(
				"CandidateSplitting: spectralIndex <= 0 !");
	}

	double Emax = Emin*pow(10, factor); 
	int n_bins = factor + 1;
	setEnergyBins(Emin, Emax, n_bins, true);
	// to compensate for loss of particles per energy bin:
	setNsplit((int) pow(10, SpectralIndex-1)); 

}


void CandidateSplitting::process(Candidate *c) const {

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
				//std::cout<< " current and previous in same energy bin " << std::endl;
				return;
			}
			//std::cout << "previous: " << prevE/GeV << " current: " << currE/GeV << std::endl;
			// current energy is in energy bin [i,i+1] or higher -> particle splitting for each crossing
			for (size_t j = i; j < Ebins.size(); ++j ){

				// adapted from Acceleration Module:
				c->updateWeight(1. / n_split); // * 1/n_split
				//std::cout<< " splitting with new weight: "  << c->getWeight() << std::endl;

				for (int i = 1; i < n_split; i++) {
				
					ref_ptr<Candidate> new_candidate = c->clone(false);
					new_candidate->parent = c;
					uint64_t snr = Candidate::getNextSerialNumber();
					new_candidate->setSerialNumber(snr);
					new_candidate->previous.setEnergy(currE); // so that new candidate is not split again in next step!
					//InteractionTag is PRIM, physically no new particles are created
					c->addSecondary(new_candidate);
					Candidate::setNextSerialNumber(snr + 1);
					//std::cout<< "new serial number" << snr << std::endl;
				}


				if (j < Ebins.size()-1 && currE < Ebins[j+1]){
					// candidate is in energy bin [j, j+1] -> no further splitting
					//std::cout<< " last bin reached, no splitting " << std::endl;
					return;
				}
			}

			//std::cout<< " return after loop " << std::endl;
			return;

		}
	}
}
	

void CandidateSplitting::setEnergyBins(double Emin, double Emax, double n_bins, bool log) {

	Ebins.resize(0);

	if (Emin > Emax){
		throw std::runtime_error(
				"CandidateSplitting: Emin > Emax!");
	}

	double dE = (Emax-Emin)/n_bins;
	//std::cout << " energy bins: " << std::endl;
	
	for (size_t i = 0; i < n_bins; ++i) {
		if (log == true) {
			Ebins.push_back(Emin * pow(Emax / Emin, i / (n_bins - 1.0)));
			//std::cout << Emin * pow(Emax / Emin, i / (n_bins - 1.0))/GeV << std::endl;
		} else {
			Ebins.push_back(Emin + i * dE);
			//std::cout << (Emin + i * dE)/GeV << std::endl;
		}
	}

}

const std::vector<double>& CandidateSplitting::getEnergyBins() const {
	return Ebins;
}

void CandidateSplitting::setNsplit(int n) {

	n_split = n;
}

int CandidateSplitting::getNsplit() const {

	return n_split;
}


} // end namespace crpropa

