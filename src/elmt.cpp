
#include "myvector.h"
#include "macros.h"
#include "elmt.h"

using namespace std;

void particle::elmtShow() const {
    cout<<"Element number "<<index<<" with "<<size<<" particles of radius "<<r<<"\n";
    cout<<"Position: "; x0.show();
    cout<<"; Velocity: "; x1.show(); cout<<";\n";
    cout<<"Orientation: "; q0.show(); cout<<";\n";
    cout<<"Inertia: " << I <<";\n";
    cout<<"Mass: "<<m<<";\n";
}

void particle::initialize(const double& partDensity, tVect& demF, unsigned int& globalIndex) {

    const double singleMass=4.0/3.0*partDensity*M_PI*r*r*r;
    m=size*singleMass;

	dens = partDensity;

    // inertia of single spheres
    I=size*2.0/5.0*singleMass*r*r;

    // translational degrees of freedom
	x1.reset();

    // rotational degrees of freedom
	w1.reset();

	wGlobal=w1;

	particleIndex=globalIndex;

    // initialize forces
	   FTotal.reset();
	   FRep.reset();
     FHydro.reset();
     FContact.reset();
     FLub.reset();
     FGrav=demF*m;
	   MTotal.reset();
     MHydro.reset();
     MContact.reset();
     MLub.reset();
	   SHydro.reset();
	   SLub.reset();

	   globalIndex++;

}

void particle::resetVelocity() {
    // resets translational and rotational velocity to zero
    // translational
    x1.reset();
    // rotational
    wGlobal.reset();
	  w1.reset();

}
