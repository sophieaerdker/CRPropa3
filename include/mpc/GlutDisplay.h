#ifndef GLUTDISPLAY_H_
#define GLUTDISPLAY_H_

#include "mpc/ModuleChain.h"

namespace mpc {

class GlutDisplay: public Module {
public:
	int counter;

	GlutDisplay();
	~GlutDisplay();

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	std::string getDescription() const;
};

} // namspace mpc

#endif /* GLUTDISPLAY_H_ */

