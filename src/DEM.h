
#ifndef DEM_H
#define	DEM_H



#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
//
#include "myvector.h"
#include "elmt.h"
#include "utils.h"
#include "LB.h"

using namespace std;

// forward declaration of IO
class IO;

class DEM{

    //IO *io;

private:
    // domain size: DEM grid is orthogonal

    // neighbor list (here COMMENTS are needed)
    // neighbor parameters
    double maxDisp;
public:
	 doubleList demSize;
	 wallList walls;
    // switchers for rotating local system
    double demInitialRepeat;
    // time step
    unsigned int demTimeStep;
  double springDamp;
  double epsilon;
  double k13;
  int nonLinType;
	int numContact;
	int numTriangle;
  int objWaitTime;
  int intervalTime;
//	int numVertex;

	double vol;

	double cPointIx[70000];
	double cPointIy[70000];
	double cPointIz[70000];
	double cPointJx[70000];
	double cPointJy[70000];
	double cPointJz[70000];
	double contactForceN[70000];


	double vertexIx[30000];
	double vertexIy[30000];
	double vertexIz[30000];

	double vertexJx[30000];
	double vertexJy[30000];
	double vertexJz[30000];

	double vertexKx[30000];
	double vertexKy[30000];
	double vertexKz[30000];

/*	int vertexI[15000];
	int vertexJ[15000];
	int vertexK[15000]; */

	int birthTime[70000];

	double avgForce;
	// stress tensor (Cauchy)

	tMat totalStress;
	tMat contactStress;
	tMat lubStress;
	tMat electroStress;
	tMat hydroStress;

//	double wallForce;

    double deltat;
    // total time duration
    double demTime;
    // force field (if not constant need to be calculated inside cycle) // should not be here!!!!!
    tVect demF;
    // rotation of the local reference system

    material sphereMat;
    // relative difference between time step in DEM and time step in LBM
    unsigned int multiStep;
    // ratio between time step and estimated duration of contacts

    particleList particles;
    // list with objects

	  objectList objects;

    bodyList bodies;

    double minPartRadius, maxPartRadius, meanPartRadius;
    // list with ghost elements

    // true if indices for the LB need to be regenerated
    bool newNeighborList;

	//contact duration

	double contactDuration;


public:
    DEM(){
        demInitialRepeat=0.0;
        demTimeStep=0;
        objWaitTime = 0;
  	  	contactDuration=0.0;
        newNeighborList=false;
        demSize.resize(3);
        demSize[0]=demSize[1]=demSize[2]=1.0;
        demTime=0.0;
        deltat=1.0;
	//	wallForce = 0.0;
        demF.reset();

        walls.clear();

        particles.clear();

        // neighbor list variables
        maxDisp=0.0;

        multiStep=1;

		totalStress.reset();
		contactStress.reset();
		lubStress.reset();
		electroStress.reset();
		hydroStress.reset();



}
    void discreteElementStep(const unsIntList& externalBoundary);
    void discreteElementGet(GetPot& config_file, GetPot& command_line);
    void discreteElementInit(const unsIntList& externalBoundary, const doubleList& externalSize, const vecList& externalBoundaryLocation, const double& externalTimeStep, const tVect externalAccel, const double& height, bool& surface);

private:
    // initialization functions

    void initializeWalls(const unsIntList& boundary, const vecList& boundaryLocation, unsigned int& globalIndex);

	// neighbor list
	void evalMaxDisp();
	void neighborList();

	// force evaluations
  void bodyObjectInteraction();
  tVect calculateEdgeDistance(tVect rb, double w2, double h2, double d2);
  void evaluateForces(int timeStep);
	void calculateParticleContactForces();
	void calculateSimplicialComplex();
	void calculateWallContactForces();
	void calculateParticleLubricationForces();
	void calculateWallLubricationForces();
	void calculateObjectParticleContactForces();
	void calculateObjectWallContactForces();
	void calculateObjectParticleLubricationForces();
	void calculateObjectWallLubricationForces();
	void calculateParticleElectroForces();
	void calculateObjectElectroForces();
	void calculateWallElectroForces();
  void calculateFictionalWallContactForces();
	// integration functions

	void integration(int timeStep);

	// periodic boundary conditions

	void xPbcs();
	void yPbcs();
	void zPbcs();

};

#endif	/* DEM_H */
