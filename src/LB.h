/*
 * File:   LB.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:06 PM
 */

#ifndef LB_H
#define LB_H

//class IO;

#include <stdio.h>
//#include <iostream>
//#include <vector>
#include <stdlib.h>
#include <math.h>
//#include <algorithm>
#include <string.h>
//
#include "macros.h"
#include "myvector.h"
#include "lattice.h"
#include "node.h"
#include "utils.h"
#include "elmt.h"
#include "getpot.h"

using namespace std;
extern ProblemName problemName;

class LB{
    // Characteristics ///////////////////////////////
    // lattice-Boltzmann equation is solved in a dimensionless space
    // with lattice spacing and time step both unitary
public: //private
	  double fluidHeight;
    // model for the fluid
    FluidMaterial fluidMaterial;
    // slip coefficient
    double slipCoefficient;
    // initial velocity for the fluid
    tVect initVelocity;
    // force field (if not constant need to be calculated inside cycle)
    tVect lbF;
    // rotation speed of the local coordinate system
    tVect lbRot;
    // cebnter of rotation of the local coordinate system
    tVect lbRotCenter;
    // total number of nodes
    unsigned int totNodes;
    // total mass (initial)
    double totalMass;
    // standard neighbors shifting
    intList ne;
    unsIntList shift;
    intList domain;
    // node types (boundary are located halfway between nodes)
    typeList types;
    //vecList positions;
    // curved info container
    //curveList curves;
    // boundary conditions
    unsIntList boundary;
    // list with nodes
    nodeList nodes;
    // indices of active nodes
    unsIntList activeNodes;
    // indices of fluid nodes
    unsIntList fluidNodes;
    // indices of solid nodes
    unsIntList particleNodes;
    // indices of interface nodes
    unsIntList interfaceNodes;
    // indices of wall nodes (only those in potential contact with the fluid)
    unsIntList wallNodes;

	unsIntList objectNodes;
public:
    // switcher for restart
    bool lbRestart;
    // restart file
     string lbRestartFile;
     // switcher for topography
    bool lbTopography;
    // topography file
     string lbTopographyFile;
     // topography container
     topography lbTop;
    // switchers for force field, non-Newtonian and everything
    bool freeSurface;
    bool forceField;
    // number of LBM step before activating the free surface
    unsigned int lbmInitialRepeat;
    // absolute time
    unsigned int time;
    // lbm size in cell units
    unsIntList lbSize;
    // lbm size in physical units (with boundaries))
    doubleList lbPhysicalSize;
    // lbm size in physical units (without boundaries))
    doubleList lbInnerPhysicalSize;
    // lbm box boundary locations
    vecList lbBoundaryLocation;
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    measureUnits unit;
    // problem-specific stuff ///////////////////////////////
    // stuff for shear cell
    double maxVisc, maxPlasticVisc, maxYieldStress, maxShearRate;
    unsigned int viscSteps, shearRateSteps;
    // stuff for drum
    double fluidMass;
    // stuff for avalanches (net-like)
    double avalanchePosit;
    // standard constructor
public:
    LB(){
        //
        freeSurface=false;
        forceField=false;
        lbRestart=false;
        //
        time=0;
        lbSize.resize(3);
        lbSize[0]=lbSize[1]=lbSize[2]=1;
        //
        slipCoefficient=0.0;
        //
        initVelocity.reset();
        lbF.reset();
        //
        totNodes=0;
        //
        shift.resize(3);
        shift[0]=shift[1]=shift[2]=1;
        domain.resize(3);
        domain[0]=domain[1]=domain[2]=1;
        //
        boundary.resize(6);
        ne.resize(lbmDirec);
        // initialize variable containers
        nodes.clear();
        types.clear();
        // initialize lists
        activeNodes.clear();
        fluidNodes.clear();
        particleNodes.clear();
        interfaceNodes.clear();
        wallNodes.clear();
		objectNodes.clear();
        //
        lbmInitialRepeat=0;
        // stuff for shear cell
        maxVisc=0.0;
        maxPlasticVisc=0.0;
        maxYieldStress=0.0;
        maxShearRate=0.0;
        viscSteps=0;
        shearRateSteps=0;
        // stuff for drum
        fluidMass=0.0;
        }
    // showing Lattice characteristics
    void LBShow() const;
    void latticeDefinition();
    void latticeBoltzmannGet(GetPot& lbmCfgFile, GetPot& command_line);
    void latticeBolzmannInit(wallList& walls, particleList& particles,objectList& objects);
    void latticeBolzmannStep(particleList& particles, wallList& walls, objectList& objects);
    void latticeBoltzmannCouplingStep(bool& newNeighborList, particleList& particles, objectList& objects);
    void latticeBoltzmannFreeSurfaceStep();

private:
    // node creation and deletion
    void createNode(const unsigned int& it);
    void deleteNode(const unsigned int& it);
    // initialize neighbors of a node coherently with the lattice
    void findNeighbors(const int& it, unsigned int* d);
    //void initializeNodeNeighbors(const int& it, node*& nodeHere);
    // initialization functions
    void initializeTypes(particleList& particles, wallList& walls , objectList& objects);
    void initializeNodes();
    void initializeLatticeBoundaries();
    void initializeParticleBoundaries(particleList& particles);
    void initializeWallBoundaries(wallList& walls);
    void initializeTopography();
    void initializeObjectBoundaries(objectList& objects);
	  void updateObjectBoundaries(objectList& objects);
    void initializeInterface();
    void restartInterface(ifstream& fluidFileID) ;
    void initializeLists();
    void resetLists();
    void initializeVariables();
    void restartVariables(ifstream& fluidFileID);
    void initializeWalls( wallList& walls);
    // integration functions
    void reconstruction(particleList& particles);
    void collision();
    void streaming(objectList& objects, wallList& walls);
    // error dumping functions
    void cleanLists();
    // coupling functions
    void updateLists(unsIntList& newPopUpNodes, unsIntList& newSolidNodes);
    void computeHydroForces(particleList& particles);
    void findNewActive(unsIntList& newPopUpNodes,  particleList& particles);
	  void findNewActiveObject(objectList& objects);
    void findNewSolid(unsIntList& newSolidNodes,  particleList& particles);
    void activeToSolid(unsIntList& newSolidNodes,  double& massSurplus);
    void solidToActive(unsIntList& newPopUpNodes, double& massSurplus);
    void updateBoundary(unsIntList& newPopUpNodes, unsIntList& newSolidNodes);
    void updateIndices(particleList& particles);
//    void updateBoundaryOld(intList& newPopUpNodes, intList& newSolidNodes, particle& dummyParticle);
    // interface functions
    void updateMass();
    void updateInterface();
    void findInterfaceMutants(unsIntList& filledNodes, unsIntList& emptiedNodes);
    void smoothenInterface(unsIntList& newFluidNodes, unsIntList& newGasNodes, unsIntList& newInterfaceNodes, double& massSurplus);
    void updateMutants(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus);
    void removeIsolated(double& massSurplus);
    void redistributeMass(const double& massSurplus);
    void enforceMassConservation();
public:
    // function for linearized index management
    tVect getPosition(const unsigned int& index) const;
    unsigned int getX(const unsigned int& index) const;
    unsigned int getY(const unsigned int& index) const;
    unsigned int getZ(const unsigned int& index) const;

};

#endif	/* LB_H */
