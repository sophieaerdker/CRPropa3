#ifndef CRPROPA_OBLIQUESHOCKMAGNETICFIELD_H
#define CRPROPA_OBLIQUESHOCKMAGNETICFIELD_H

#include "crpropa/magneticField/MagneticField.h"


#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "crpropa/Vector3.h"
#include "crpropa/Referenced.h"
#include "crpropa/Units.h"

namespace crpropa {

/**

@class Oblique Shock Magnetic Field
@brief Magnetic field model for oblique shocks, use together with ObliqueAdvectionShock!

*/

class ObliqueShockMagneticField: public MagneticField {
private:
	double Bx_up; // upstreamn B -component perp to shock front
	double By; // constant B-component parallel to shock front
	double r_comp; // shock compression ratio
	double x_sh; // shock widht

public:
/** Constructor
	@param  Bx_up // upstreamn B -component perp to shock front
	@param By // constant B-component parallel to shock front
	@param r_comp // shock compression ratio
	@param x_sh // shock widht
*/
	ObliqueShockMagneticField(double Bx_up, double By, double r_comp, double xsh);

	Vector3d getField(const Vector3d &pos) const;	
		
	void setBx(double Bx_up);
	void setBy(double By);
	void setComp(double r_comp);
	void setShockwidth(double x_sh);

};

	 
} // end namespace crpropa

#endif // CRPROPA_ACHIMEDEANSPIRALFIELD_H
