

#ifndef ELMT_H
#define	ELMT_H

#include "myvector.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/




class particle{
public:
    unsigned int index;
    // constitutive particles indexes

	  unsigned int particleIndex;

    unsigned int size;
    // radius of constitutive particles (supposed, for the moment, of being constant)
    double r;
    // mass
    double m;
    // moment of inertia
    double I;

    double coord;
    double lub_coord;

	  double dens;
    // Position and derivatives
    // position of the center of mass of the element
    tVect  x0;
    tVect x_i;

    // velocity of the center of the element
    tVect x1;

	//neighbor list parameters

	int    list[50];
	int    _list[50];
	double xix[50];
	double xiy[50];
	double xiz[50];
	double _xix[50];
	double _xiy[50];
	double _xiz[50];

	tVect xiw;
	tVect xio;
	tVect x0_lst;

    // for verlet algorithm
    // mass
    double _m;
    // moment of inertia
    double _I;
    tVect _x0;
	tVect xt0;
	tVect disp;

	double disp_r;

	double trueStrain;
	double engStrain;
	double trueMod;
	double engMod;

	double mu;
	double ratio;

	double pressure;
	double normal1;
	double normal2;

    // Orientation and derivatives
    // orientation of the center of mass of the element
    tQuat q0, q1;
    // quaternion rates (related to angular velocity, acceleration...)
    // angular velocity rates (global)

    tVect w0,w1,wGlobal;

    // forces and moments
    tVect FHydro,FContact,FTotal,FLub,FGrav,FRep;
    tVect FElas, FVisc, FVisc_h;
    tVect MHydro,MContact,MTotal,MLub;
	  tMat SHydro, SLub, SContact, SRep, STotal;
    tMat SElas, SVisc, SVisc_h;
    tVect overlap;

    // fluid mass entrained by the particle
    double fluidVolume;


	particle() {

		index=0;
		particleIndex =  0;
        size=1;
        r=1.0;
        m=1.0;
        I=1.0;
		dens=1.0;

        disp=x0=x1=xiw=xio=x0_lst=xt0=_x0=tVect(0.0,0.0,0.0);

		pressure = normal1 = normal2 = disp_r = 0.0;


		q0=tQuat(1.0,0.0,0.0,0.0);
        q1=tQuat(0.0,0.0,0.0,0.0);

        w0=w1=tVect(0.0,0.0,0.0);

		list[0] = 0;
		_list[0]= 0;
		xix[0] = 0.0;
		xiy[0] = 0.0;
		xiz[0] = 0.0;
		_xix[0] = 0.0;
		_xiy[0] = 0.0;
		_xiz[0] = 0.0;

		FGrav.reset();
		FHydro.reset();
    FContact.reset();
    FLub.reset();
    FTotal.reset();
		FRep.reset();
    FElas.reset();
    FVisc.reset();

		SHydro.reset();
		SLub.reset();
		SContact.reset();
		SRep.reset();
		STotal.reset();
    SElas.reset();
    SVisc.reset();


		MHydro.reset();
    MContact.reset();
    MLub.reset();
    MTotal.reset();

    fluidVolume=0.0;

    }
    void elmtShow()const;
    void initialize(const double& partDensity, tVect& demF, unsigned int& globalIndex);
    void resetVelocity();
};

class object{
public:
    // object index
    unsigned int index;
    // submerged volume
    double subVol;

    double wPlate;
    double hPlate;
    double dPlate;
    double r;
    //double omega;
    //deepest point of the impactor
    tVect dPoint;
	  double m;
	  double _m;
    double rho;
	  double I;
	  double _I;
    tMat Ibody;
    tMat _Ibody;
    tMat Rot;
    tMat RotT;
    tQuat w0;
    tQuat _w0;
    // position of the object
    tVect x0;
	  tVect _x0;
    // position of each vertex of the box
    tVect p1;
    tVect p2;
    tVect p3;
    tVect p4;
    tVect p5;
    tVect p6;
    tVect p7;
    tVect p8;
    // velocity of the object
    tVect x1, w1;
    // force on the object
    tVect FContact,FHydro, FLub, FGrav, FRep, FSpring, FTotal;
	  tVect MContact, MHydro, MTotal, MLub;
    tVect FElas, FVisc, FVisc_h;

    object() {
        index=0;
        wPlate=0.0;
        hPlate=0.0;
        dPlate=0.0;
        subVol = 0.0;
		    _m=m=1.0;
        _x0=x0=tVect(0.0,0.0,0.0);
        x1=tVect(0.0,0.0,0.0);
        FContact=tVect(0.0,0.0,0.0);
        FSpring=tVect(0.0,0.0,0.0);
        FHydro=tVect(0.0,0.0,0.0);
		    FGrav=tVect(0.0,0.0,0.0);
		    FRep=tVect(0.0,0.0,0.0);
		    FLub=tVect(0.0,0.0,0.0);
        FElas=tVect(0.0,0.0,0.0);
        FVisc=tVect(0.0,0.0,0.0);
		    FTotal=tVect(0.0,0.0,0.0);
		    MTotal=tVect(0.0,0.0,0.0);
		    MContact=tVect(0.0,0.0,0.0);
	      MLub=tVect(0.0,0.0,0.0);
		    MHydro=tVect(0.0,0.0,0.0);
        _w0 = w0 = tQuat(1.0,0.0,0.0,0.0);
    }
};


class body{
public:
    // object index
    unsigned int index;
    double r;
	  double m;
    double _m;
    // position of the object
    tVect x0;
    double L0;
    double L;
    double k;
    tVect x1;
    // force on the object
    tVect FGrav, FSpring, FTotal;
    body() {
        index=0;
		    m=1.0;
        r=1.0;
        L=L0=1.0;
        k=1.0;
        x0=tVect(0.0,0.0,0.0);
        x1=tVect(0.0,0.0,0.0);
        FSpring=tVect(0.0,0.0,0.0);
		    FGrav=tVect(0.0,0.0,0.0);
		    FTotal=tVect(0.0,0.0,0.0);
    }
};


class material{
public:
    // density [ mass / length² ]
    double density;

	  double objectDensity;

    double bodyDensity;

    // linear model ////////////////////
    // stiffness
    double linearStiff;

    // normal damping ///////////////////////////
    // normal viscous coefficient
    double dampCoeff;

    // tangential model ////////////////////////
    // particle particle friction
    double frictionCoefPart;
    // particle-wall friction

	double lubNormal;
	double lubTangential;
	double lubRotation;

	double lubCutoff;
	double kinVisc;

	double electroForce;
	double debyeLength;


    // default constructor
    material(){
        density=1.0;
        linearStiff=0.0;
        dampCoeff=0.0;
        frictionCoefPart=0.0;
		lubNormal = 0.0;
		lubTangential = 0.0;
		lubRotation = 0.0;
		lubCutoff = 0.0;
		kinVisc = 0.0;
    }
};



#endif	/* ELMT_H */
