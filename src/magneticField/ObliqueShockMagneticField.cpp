#include "crpropa/magneticField/ObliqueShockMagneticField.h"

namespace crpropa {

ObliqueShockMagneticField::ObliqueShockMagneticField(double Bx_up, double By, double r_comp, double xsh) {

	// magnetic field that follows ObliqueAdvectionShock
	setBx(Bx_up);
	setBy(By);
	setComp(r_comp);
	setShockwidth(x_sh);
}
   
   
Vector3d ObliqueShockMagneticField::getField(const Vector3d &pos) const {
	
	double x = pos.x;
	double Bx_down = Bx_up/r_comp;

	double a = (Bx_up + Bx_down)*0.5;
	double b = (Bx_up - Bx_down)*0.5;

	Vector3d B(0.);

	B.x = a - b*tanh(x/x_sh);
	B.y = By;

	return B;
}

void ObliqueShockMagneticField::setBx(double b) {
	Bx_up = b;
	return;
}

void ObliqueShockMagneticField::setBy(double b) {
	By = b;
	return;
}

void ObliqueShockMagneticField::setComp(double r) {
	r_comp = r;
	return;
}

void ObliqueShockMagneticField::setShockwidth(double w) {
	x_sh = w;
	return;
}

} //end namespace crpropa
