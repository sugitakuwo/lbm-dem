
#include "LB.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PUBLIC FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

void LB::LBShow() const {
    cout << "LATTICE CHARACTERISTICS (lattice units):" << endl;
    cout << "directions=" << lbmDirec << ";" << endl;
    cout << "time unit = " << unit.Time << "; length unit = " << unit.Length << "; density unit = " << unit.Density << ";" << endl;
    cout << "dynVisc unit = " << unit.DynVisc << "; kinVisc unit = " << unit.KinVisc << "; speed unit = " << unit.Speed << ";" << endl;
    cout << "xdim =" << lbSize[0] << "; ydim= " << lbSize[1] << "; zdim= " << lbSize[2] << ";" << endl;
    cout << "Init tau = " << (0.5 + 3.0 * fluidMaterial.initDynVisc) << ";init visc =" << fluidMaterial.initDynVisc << "; init density = " << fluidMaterial.initDensity << ";" << endl;
    cout << "Min tau = " << minTau << "; Max tau = " << maxTau << ";" << endl;
    switch (fluidMaterial.rheologyModel) {
        case BINGHAM:
        {
            cout << "Plastic viscosity = " << fluidMaterial.plasticVisc << "(tau = " << (0.5 + 3.0 * fluidMaterial.plasticVisc) << " ); yield stress=" << fluidMaterial.yieldStress << ";" << endl;
            break;
        }
        case FRICTIONAL:
        {
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            break;
        }
        case VOELLMY:
        {
            cout << "Voellmy fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << "; Collisional/turbulent constant csi=" << fluidMaterial.voellmyCoefficient << ";" << endl;
            break;
        }
    }
    cout << "Initial Velocity=";
    initVelocity.show();
    cout << ";" << endl;
    cout << "F=";
    lbF.show();
    cout << endl;
    cout << "Rotation speed =";
    lbRot.show();
    cout << ";" << endl;
    cout << "Rotation center =";
    lbRotCenter.show();
    cout << ";" << endl;
    cout << "X(-)=" << boundary[0] << "; X(+)=" << boundary[1] << "; Y(-)=" << boundary[2] << "; Y(+)=" << boundary[3] << "; Z(-)=" << boundary[4] << "; Z(+)=" << boundary[5] << ";" << endl;
    cout << "PHYSICAL CHARACTERISTICS (physical units)" << endl;
    cout << "xdim=" << (lbSize[0] - 2) * unit.Length << "; ydim= " << (lbSize[1] - 2) * unit.Length << "; zdim= " << (lbSize[2] - 2) * unit.Length << ";" << endl;
    cout << "init kin viscosity=" << fluidMaterial.initDynVisc * unit.KinVisc << ", init dyn viscosity=" << fluidMaterial.initDynVisc * unit.DynVisc << "; init density=" << fluidMaterial.initDensity * unit.Density << ";" << endl;
    switch (fluidMaterial.rheologyModel) {
        case BINGHAM:
        {
            cout << "plastic kin viscosity=" << fluidMaterial.plasticVisc * unit.KinVisc << ", plastic dyn viscosity=" << fluidMaterial.plasticVisc * unit.DynVisc << "; yield stress=" << fluidMaterial.yieldStress * unit.Stress << ";" << endl;
            break;
        }
        case FRICTIONAL:
        {
            cout << "Fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << ";" << endl;
            break;
        }
        case VOELLMY:
        {
            cout << "Voellmy fluid friction coefficient = " << fluidMaterial.frictionCoefFluid << "; Collisional/turbulent constant csi=" << fluidMaterial.voellmyCoefficient * unit.Time << ";" << endl;
            break;
        }
    }
    cout << "Min dyn viscosity = " << (minTau - 0.5) / 3.0 * unit.DynVisc << "; Max dym viscosity = " << (maxTau - 0.5) / 3.0 * unit.DynVisc << ";" << endl;
    cout << "Initial Velocity=";
    const tVect initVelocityPhysic = initVelocity * unit.Speed;
    initVelocityPhysic.show();
    cout << ";" << endl;
    cout << "F=";
    const tVect FPhysic = lbF * unit.Accel;
    FPhysic.show();
    cout << ";" << endl;
    cout << "Rotation speed =";
    const tVect RotPhysic = lbRot * unit.AngVel;
    RotPhysic.show();
    cout << ";" << endl;
    cout << "Rotation center =";
    const tVect RotCenterPhysic = lbRotCenter * unit.Length;
    RotCenterPhysic.show();
    cout << ";" << endl;
}

void LB::latticeDefinition() {
    // LATTICE PARAMETERS  ////////////////
    //size-dependent; the others are in the lattice.h

    // domain movement variables
    shift[2] = lbSize[0] * lbSize[1];
    shift[1] = lbSize[0];
    shift[0] = 1;
    //
    domain[2] = lbSize[2] * lbSize[1] * lbSize[0] - 2 * shift[2];
    domain[1] = lbSize[1] * lbSize[0] - 2 * shift[1];
    domain[0] = lbSize[0] - 2 * shift[0];

    // standard neighbors shifting
    // order is O,x,y,z,xy,yz,zx.
    ne[0] = 0;
    //
    ne[1] = shift[0];
    ne[2] = -shift[0];
    //
    ne[3] = shift[1];
    ne[4] = -shift[1];
    //
    ne[5] = shift[2];
    ne[6] = -shift[2];
    //
    ne[7] = shift[0] + shift[1];
    ne[8] = -shift[0] - shift[1];
    ne[9] = -shift[0] + shift[1];
    ne[10] = shift[0] - shift[1];
    //
    ne[11] = shift[1] + shift[2];
    ne[12] = -shift[1] - shift[2];
    ne[13] = -shift[1] + shift[2];
    ne[14] = shift[1] - shift[2];
    //
    ne[15] = shift[2] + shift[0];
    ne[16] = -shift[2] - shift[0];
    ne[17] = -shift[2] + shift[0];
    ne[18] = shift[2] - shift[0];

}

void LB::latticeBoltzmannGet(GetPot& configFile, GetPot& commandLine) {


    // conversion units //////////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // measure units for DEM solver are in the international system

    cout << "Getting LBM info from config file" << endl;

    PARSE_CLASS_MEMBER(configFile, lbRestart, "lbRestart", false);
    PARSE_CLASS_MEMBER(configFile, lbRestartFile, "lbRestartFile", "");
    PARSE_CLASS_MEMBER(configFile, lbTopography, "lbTopography", false);
    PARSE_CLASS_MEMBER(configFile, lbTopographyFile, "lbTopographyFile", "");
    if (lbTopography && !lbRestart) {
        cout << "Using topography file for boundary and initial surface:  " << lbTopographyFile << endl;
    } else if (lbTopography && lbRestart) {
        cout << "Using topography file for boundary: " << lbTopographyFile << ", and restart file for fluid:  " << lbRestartFile << endl;
    } else if (lbRestart) {
        cout << "Computation from LB restart file " << lbRestartFile << endl;
    } else {
        cout << "Computation from scratch - for interface geometry check LB::initializeInterface() " << endl;
    }

    // primary
    PARSE_CLASS_MEMBER(configFile, unit.Length, "unitLength", 1.0);
    ASSERT(unit.Length > 0);
    PARSE_CLASS_MEMBER(configFile, unit.Time, "unitTime", 1.0);
    //ASSERT(unit.Time > 0);
    PARSE_CLASS_MEMBER(configFile, unit.Density, "unitDensity", 1.0);
    ASSERT(unit.Density > 0);

    // fixing length: domain size
    double lbSizeX, lbSizeY, lbSizeZ;
    lbSize.resize(lbmDim);
    PARSE_CLASS_MEMBER(configFile, lbSizeX, "lbSizeX", 0.0);
    ASSERT(lbSizeX > 0.0);
    PARSE_CLASS_MEMBER(configFile, lbSizeY, "lbSizeY", 0.0);
    ASSERT(lbSizeY > 0.0);
    PARSE_CLASS_MEMBER(configFile, lbSizeZ, "lbSizeZ", 0.0);
    ASSERT(lbSizeZ > 0.0);
    // scaling and adapting length unit to get exact domain discretization
    lbSize[0] = int(floor(lbSizeX / unit.Length + 0.5));
    lbSize[1] = int(floor(lbSizeY / unit.Length + 0.5));
    lbSize[2] = int(floor(lbSizeZ / unit.Length + 0.5));

    // domain in physical units
    lbPhysicalSize.resize(3);
    lbPhysicalSize[0] = double(lbSize[0]) * unit.Length;
    lbPhysicalSize[1] = double(lbSize[1]) * unit.Length;
    lbPhysicalSize[2] = double(lbSize[2]) * unit.Length;

    // domain inside boundaries
    lbInnerPhysicalSize.resize(3);
    lbInnerPhysicalSize[0] = double(lbSize[0] - 2.0) * unit.Length;
    lbInnerPhysicalSize[1] = double(lbSize[1] - 2.0) * unit.Length;
    lbInnerPhysicalSize[2] = double(lbSize[2] - 2.0) * unit.Length;

    // location of the domain boundary
    lbBoundaryLocation.resize(6);
    lbBoundaryLocation[0] = tVect(0.5 * unit.Length, 0.0, 0.0);
    lbBoundaryLocation[1] = tVect(double(lbSize[0] - 1.5) * unit.Length, 0.0, 0.0);
    lbBoundaryLocation[2] = tVect(0.0, 0.5 * unit.Length, 0.0);
    lbBoundaryLocation[3] = tVect(0.0, double(lbSize[1] - 1.5) * unit.Length, 0.0);
    lbBoundaryLocation[4] = tVect(0.0, 0.0, 0.5 * unit.Length);
    lbBoundaryLocation[5] = tVect(0.0, 0.0, double(lbSize[2] - 1.5) * unit.Length);


    // material //////////////////////////////////////////////////
    string rheologyModelString;
    PARSE_CLASS_MEMBER(configFile, rheologyModelString, "rheologyModel", "none");
    if (rheologyModelString == "NEWTONIAN") fluidMaterial.rheologyModel = NEWTONIAN;
    else if (rheologyModelString == "BINGHAM") fluidMaterial.rheologyModel = BINGHAM;
    else if (rheologyModelString == "FRICTIONAL") fluidMaterial.rheologyModel = FRICTIONAL;
    else if (rheologyModelString == "VOELLMY") fluidMaterial.rheologyModel = VOELLMY;

    switch (fluidMaterial.rheologyModel) {
        case NEWTONIAN:
        {
            // getting viscosity
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.initDynVisc, "initVisc", 1.0);
            ASSERT(fluidMaterial.initDynVisc > 0.0);
            break;
        }
        case BINGHAM:
        {
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.plasticVisc, "plasticVisc", 0.0);
            ASSERT(fluidMaterial.plasticVisc >= 0.0);
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.yieldStress, "yieldStress", 0.0);
            ASSERT(fluidMaterial.yieldStress >= 0.0);
            break;
        }
        case FRICTIONAL:
        {
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.frictionCoefFluid, "frictionCoefFluid", 0.0);
            ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
            // getting initial viscosity
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.initDynVisc, "initVisc", 1.0);
            ASSERT(fluidMaterial.initDynVisc > 0.0);
            break;
        }
        case VOELLMY:
        {
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.frictionCoefFluid, "frictionCoefFluid", 0.0);
            ASSERT(fluidMaterial.frictionCoefFluid >= 0.0 && fluidMaterial.frictionCoefFluid <= 1.0);
            // getting initial viscosity
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.voellmyCoefficient, "voellmyConstant", 1.0);
            ASSERT(fluidMaterial.voellmyCoefficient >= 0.0);
            // getting initial viscosity
            PARSE_CLASS_MEMBER(configFile, fluidMaterial.initDynVisc, "initVisc", 1.0);
            ASSERT(fluidMaterial.initDynVisc > 0.0);
            break;
        }
    }
    // are we using turbulence modeling? (shutting down improves performances)
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.turbulenceOn, "turbulenceSolve", 0);
    PARSE_CLASS_MEMBER(configFile, fluidMaterial.turbConst, "turbConst", 0.0);
    ASSERT(fluidMaterial.turbConst >= 0.0);

    // if time step is chosen as optimal, compute it
    // see thesis, §4.1.1
    if (unit.Time == 0.0) {
        const double tauOpt = 1.0;
        switch (fluidMaterial.rheologyModel) {
            case BINGHAM:
            {
                unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.plasticVisc * (1.0 / 3.0)*(tauOpt - 0.5);
                break;
            }
            default:
            {
                unit.Time = unit.Density * unit.Length * unit.Length / fluidMaterial.initDynVisc * (1.0 / 3.0)*(tauOpt - 0.5);
                break;
            }
                cout << "Time step automatically chosen: " << unit.Time << " s" << endl;
        }
    }


    // compute all non-primary conversion units
    unit.setComposite();

        double initVelocityX, initVelocityY, initVelocityZ;
    PARSE_CLASS_MEMBER(configFile, initVelocityX, "initVelocityX", 0.0);
    PARSE_CLASS_MEMBER(configFile, initVelocityY, "initVelocityY", 0.0);
    PARSE_CLASS_MEMBER(configFile, initVelocityZ, "initVelocityZ", 0.0);
    initVelocity = tVect(initVelocityX, initVelocityY, initVelocityZ);
    initVelocity /= unit.Speed;

    // external force field (as acceleration)
    double lbFX, lbFY, lbFZ;
    PARSE_CLASS_MEMBER(configFile, lbFX, "lbFX", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbFY, "lbFY", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbFZ, "lbFZ", 0.0);
    lbF = tVect(lbFX, lbFY, lbFZ);


    lbF /= unit.Accel;

        // rotation of the local coordinate system (rad/s)
    double lbRotX, lbRotY, lbRotZ;
    PARSE_CLASS_MEMBER(configFile, lbRotX, "lbRotX", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbRotY, "lbRotY", 0.0);
    PARSE_CLASS_MEMBER(configFile, lbRotZ, "lbRotZ", 0.0);
    lbRot = tVect(lbRotX, lbRotY, lbRotZ);
    lbRot /= unit.AngVel;

    // position of the rotation center of the local coordinate system
    // default is (-1.0,-1.0,-1.0) to avoid computing a zero distance from the center
    double lbRotCenterX, lbRotCenterY, lbRotCenterZ;
    PARSE_CLASS_MEMBER(configFile, lbRotCenterX, "lbRotCenterX", -1.0);
    PARSE_CLASS_MEMBER(configFile, lbRotCenterY, "lbRotCenterY", -1.0);
    PARSE_CLASS_MEMBER(configFile, lbRotCenterZ, "lbRotCenterZ", -1.0);
    lbRotCenter = tVect(lbRotCenterX, lbRotCenterY, lbRotCenterZ);
    lbRotCenter /= unit.Length;


    boundary.resize(6);
    PARSE_CLASS_MEMBER(configFile, boundary[0], "boundary0", 1);
    PARSE_CLASS_MEMBER(configFile, boundary[1], "boundary1", 1);
    PARSE_CLASS_MEMBER(configFile, boundary[2], "boundary2", 1);
    PARSE_CLASS_MEMBER(configFile, boundary[3], "boundary3", 1);
    PARSE_CLASS_MEMBER(configFile, boundary[4], "boundary4", 1);
    PARSE_CLASS_MEMBER(configFile, boundary[5], "boundary5", 1);

    PARSE_CLASS_MEMBER(configFile, slipCoefficient, "slipCoefficient", 0.0);
    ASSERT(slipCoefficient >= 0.0);

    PARSE_CLASS_MEMBER(configFile, fluidHeight, "fluidHeight", 0.0);
    ASSERT(fluidHeight >= 0.0);

    // scaling rheological parameters
    fluidMaterial.initDynVisc /= unit.DynVisc;
    fluidMaterial.plasticVisc /= unit.DynVisc;
    fluidMaterial.yieldStress /= unit.Stress;
    fluidMaterial.voellmyCoefficient /= unit.Accel;
    // scaling Voellmy parameter by "gravity"
    const double voellmyGravity=lbF.norm();
    fluidMaterial.voellmyConstant=fluidMaterial.voellmyCoefficient/voellmyGravity;
    //    PARSE_CLASS_MEMBER(lbmCfgFile, initDensity, "initDensity",1.0);
    //    initDensity/=unit.Density;
    // to avoid errors we set this to be just 1
    fluidMaterial.initDensity = 1.0;
    ASSERT(fluidMaterial.initDensity > 0.0);






}

void LB::latticeBolzmannInit(wallList& walls, particleList& particles, objectList& objects) {
    //  Lattice Boltzmann initialization steps

    // first comes the initialization of the data structures
    initializeNodes();

    // application of lattice boundaries
    initializeLatticeBoundaries();

    // then the initial node type must be identified for every node (if not specified, it is already Fluid)
    initializeTypes(particles, walls, objects);

    ifstream fluidFileID;
    if (lbRestart) {
        // open fluid restart file
        fluidFileID.open(lbRestartFile.c_str(), ios::in);
        ASSERT(fluidFileID.is_open());
        // check if the restart file size is ok
        unsigned int restartX, restartY, restartZ;
        fluidFileID>>restartX;
        fluidFileID>>restartY;
        fluidFileID>>restartZ;
        ASSERT(restartX == lbSize[0]);
        ASSERT(restartY == lbSize[1]);
        ASSERT(restartZ == lbSize[2]);
        restartInterface(fluidFileID);
        // initialize variables for active nodes
        restartVariables(fluidFileID);
    } else {

	 if	(freeSurface){
        // initialize interface
        initializeInterface();
	}
        // initialize variables for active nodes
        initializeVariables();
    }

    // initialize variables for wall nodes
    initializeWalls(walls);

    // initialize node properties and lists
    initializeLists();

    //create active list and checks for errors in the list
    cleanLists();

    // in case mass needs to be kept constant, compute it here
    totalMass = 0.0;

    for (int it = 0; it < activeNodes.size(); ++it) {
        const unsigned int index = activeNodes[it];
        if (!types[index].isInsideParticle()) {
            totalMass += nodes[index]->mass;
        }
    }

    cout << "Done with initialization" << endl;
}

void LB::latticeBolzmannStep(particleList& particles, wallList& walls, objectList& objects) {

    // Lattice Boltzmann core steps

	initializeLists();

	cleanLists();

    // reconstruction of macroscopic variables
    reconstruction(particles);
    // compute interaction forces
    computeHydroForces(particles);
    //collision operator
    collision();
    //streaming operator
    streaming(objects, walls);


}



void LB::latticeBoltzmannFreeSurfaceStep() {
    // in case mass needs to be kept constant, call enforcing function here

    // mass and free surface update
    updateMass();
    updateInterface();
    cleanLists();

}

void LB::latticeBoltzmannCouplingStep(bool& newNeighborList, particleList& particles, objectList& objects) {
    // identifies which nodes need to have an update due to particle movement
    // the complexity arises from trying to keep the scaling as close to linear as possible
    // maybe the best idea is to do this in two step:
    // 1) the first is to check for new active nodes and initialize them
    // 2) the second is checking for new solid nodes.
    // this automatically check also for the hideous case of particle to particle double transition

    // first we check if a new neighbor table has been defined. In that case, the indexing needs to be reinitialized
	updateObjectBoundaries(objects);

	 if (newNeighborList) {
  	   updateIndices(particles);
	    newNeighborList = false;
   	}
    // declaring list for mutant nodes
    unsIntList newPopUpNodes, newSolidNodes;
    // emptying list of mutant nodes
    newPopUpNodes.clear();
    newSolidNodes.clear();
    // double massSurplus=0.0;

    // SOLID TO ACTIVE CHECK
    findNewActive(newPopUpNodes,  particles);

    //    solidToActive(newPopUpNodes, elmts, massSurplus);

    // ACTIVE TO SOLID CHECK
    findNewSolid(newSolidNodes, particles);

    //    activeToSolid(newSolidNodes, elmts, massSurplus);

    // redistributing extra mass due to popping of nodes
    // redistributeMass(massSurplus);

}


/*//////////////////////////////////////////////////////////////////////////////////////////////////////
 // PRIVATE FUNCTIONS
 //////////////////////////////////////////////////////////////////////////////////////////////////////*/

// node creation and deletion

void LB::createNode(const unsigned int& it) {
    // creating node
    nodes[it] = new node;
    // initialize its neighbors
    findNeighbors(it,nodes[it]->d);
    //initializeNodeNeighbors(it, nodes[it]);

}

// initialization functions

void LB::initializeTypes(particleList& particles, wallList& walls , objectList& objects) {

    // application of particle initial position
    initializeParticleBoundaries(particles);
    // application of solid walls
    initializeWallBoundaries(walls);

	initializeObjectBoundaries(objects);
    // initializing topography if one is present
    initializeTopography();

}

void LB::initializeNodes() {

    cout << "Initializing nodes containers and types" << endl;

    // total number of nodes
    totNodes = lbSize[0] * lbSize[1] * lbSize[2];

    const int maxSize=nodes.max_size();
    ASSERT(totNodes<maxSize);

    // vector with node pointers
    //nodes.reserve(totNodes);
    //nodes.resize(totNodes,0);

    nodes.resize(totNodes,0);
    #pragma omp parallel for
    for (int it = 0; it < totNodes; ++it) {
        nodes[it] = 0;
    }
    // vector with types
    //types.reserve(totNodes);
    types.resize(totNodes);
    #pragma omp parallel for
    for (int it = 0; it < totNodes; ++it) {
        types[it].setFluid();
    }
    // vector with positions
    //positions.resize(totNodes);
    //for (int it = 0; it < totNodes; ++it) {
    //    positions[it] = setPosition(it);
    //}
    // vector with curves
    //curves.resize(totNodes);
    //for (int it = 0; it < totNodes; ++it) {
    //    curves[it] = 0;
    //}
}

void LB::initializeLatticeBoundaries() {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    // BOUNDARY CONDITIONS ///////////////////////////
    // solid boundary wins over all in corners, where more than 1 bc is defined
    cout << "Initializing boundaries" << endl;
#pragma omp parallel for
    for (int it = 0; it < totNodes; ++it) {
        if (!types[it].isWall()) {
            if (getX(it) == 0) {
                types[it].setType(boundary[0]);
            } else if (getX(it) == lbSize[0] - 1) {
                types[it].setType(boundary[1]);
            } else if (getY(it) == 0) {
                types[it].setType(boundary[2]);
            } else if (getY(it) == lbSize[1] - 1) {
                types[it].setType(boundary[3]);
            } else if (getZ(it) == 0) {
                types[it].setType(boundary[4]);
            } else if (getZ(it) == lbSize[2] - 1) {
                types[it].setType(boundary[5]);
            }
        }
    }

}


void LB::findNeighbors(const int& it, unsigned int* d) {
    // assign boundary characteristic to nodes (see class)
    // if not differently defined, type is 0 (fluid)

    // resetting neighbors and vectors
    for (int j = 0; j < lbmDirec; ++j) {
        d[j] = it + ne[j];
    }

    // BOUNDARY CONDITIONS ///////////////////////////
    // nodes on the boundary have no neighbors
    if (types[it].isWall()) {
        const unsigned int xHere = getX(it);
        if (xHere == 0)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) < 0.0) {
                    d[j] = it;
                }
            }
        if (xHere == lbSize[0] - 1)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) > 0.0) {
                    d[j] = it;
                }
            }
        const unsigned int yHere = getY(it);
        if (yHere == 0)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) < 0.0) {
                    d[j] = it;
                }
            }
        if (yHere == lbSize[1] - 1)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) > 0.0) {
                    d[j] = it;
                }
            }
        const unsigned int zHere = getZ(it);
        if (zHere == 0)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) < 0.0) {
                    d[j] = it;
                }
            }
        if (zHere == lbSize[2] - 1)
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) > 0.0) {
                    d[j] = it;
                }
            }
    }
// PERIODICITY ////////////////////////////////////////////////
    // assigning periodicity conditions (this needs to be done after applying boundary conditions)
    // runs through free cells and identifies neighboring cells. If neighbor cell is
    // a special cell (periodic) then the proper neighboring condition is applied
    // calculates the effect of periodicity
    else if (types[it].isActive()) {
        // neighboring and periodicity vector for boundary update
        int pbc[lbmDirec];
        for (int j = 1; j < lbmDirec; ++j) {
            pbc[j] = 0;
        }


        if (types[d[1]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) > 0.0) {
                    pbc[j] -= domain[0];
                }
            }
        }
        if (types[d[2]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Xp) < 0.0) {
                    pbc[j] += domain[0];
                }
            }
        }
        if (types[d[3]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) > 0.0) {
                    pbc[j] -= domain[1];
                }
            }
        }
        if (types[d[4]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Yp) < 0.0) {
                    pbc[j] += domain[1];
                }
            }
        }
        if (types[d[5]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) > 0.0) {
                    pbc[j] -= domain[2];
                }
            }
        }
        if (types[d[6]].isPeriodic()) {
            for (int j = 1; j < lbmDirec; ++j) {
                if (v[j].dot(Zp) < 0.0) {
                    pbc[j] += domain[2];
                }
            }
        }

        // apply periodicity
        for (int j = 1; j < lbmDirec; ++j) {
            d[j] += pbc[j];
        }
    }
}


void LB::initializeParticleBoundaries(particleList& particles) {

    // INITIAL PARTICLE POSITION ////////////////////////
  cout << "Initializing particle nodes" << endl;

    // set all outside
#pragma omp parallel for
    for (int it = 0; it < totNodes; ++it) {
        // reset old list
        types[it].setOutsideParticle();
    }

    for (int n = 0; n < particles.size(); ++n) {
        const tVect convertedPosition = particles[n].x0 / unit.Length;
        const double convertedRadius = particles[n].r / unit.Length;
        const unsigned int indexHere = particles[n].particleIndex;
#pragma omp parallel for
        for (int ia = 0; ia < activeNodes.size(); ++ia) {
        const unsigned int index=activeNodes[ia];
            // checking if we are not in a boundary
            // whose properties have already been assigned and must be invariant in time
                // checking if node is inside a particle
                if (getPosition(index).insideSphere(convertedPosition, convertedRadius)) { //-0.5?
                    types[index].setInsideParticle();
                    types[index].setSolidIndex(indexHere);
                }
        }
    }
}

void LB::initializeWallBoundaries(wallList& walls) {

    // SOLID WALLS ////////////////////////
 //   cout << "Initializing solid walls" << endl;
    for (int iw = 0; iw < walls.size(); ++iw) {
        const tVect convertedWallp = walls[iw].p / unit.Length;
        const tVect normHere = walls[iw].n;
        const unsigned int indexHere = walls[iw].index;
        const bool slipHere = walls[iw].slip;
        const bool movingHere = walls[iw].moving;
#pragma omp parallel for
        for (int it = 0; it < totNodes; ++it) {
            // creating solid cells
            if (getPosition(it).insidePlane(convertedWallp, normHere)) {
                // setting solidIndex
                types[it].setSolidIndex(indexHere);
                // setting type: 5-6=slip, 7-8=no-slip
                if (slipHere) {
                    // setting type for slip: 5=static, 6=moving
                    if (movingHere) {
                        types[it].setSlipDynWall();
                    } else {
                        types[it].setSlipStatWall();
                    }
                } else {
                    // setting type for no-slip: 7=static, 8=moving
                    if (movingHere) {
                        types[it].setDynWall();
                    } else {
                        types[it].setStatWall();
                    }
                }
            }
        }
    }
}


void LB::initializeObjectBoundaries(objectList& objects) {
	  cout << "Initializing objects" << endl;
	const double zero = 0.0;
  tVect Zero = tVect(0.0,0.0,0.0);
	objectNodes.clear();

    for (int io = 0; io < objects.size(); io++) {
        const tVect convertedPosition = objects[io].x0 / unit.Length;
        const double convertedRadius = objects[io].r / unit.Length;
        const unsigned int indexHere = objects[io].index;
#pragma omp parallel for
        for (int it = 0; it < totNodes; ++it) {
            const tVect nodePosition = getPosition(it);
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) {// setting solidIndex
                types[it].setSolidIndex(indexHere);
                types[it].setObjectWall();
            }
        }
    }

	for (int ip = 0; ip < totNodes; ip++ ){

		if ( types[ip].isObjectWall()) {

			createNode(ip);

        	tVect solidVelocity;
        	const tVect nodePosition = getPosition(ip);
        	unsigned int solidIndex = types[ip].getSolidIndex();

        	if (nodePosition.insideSphere(objects[solidIndex].x0 / unit.Length, objects[solidIndex].r / unit.Length)) {
            	solidVelocity = objects[solidIndex].x1 / unit.Speed;

        // reset velocity and mass (useful for plotting)
        // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
        		nodes[ip]->initialize(fluidMaterial.initDensity, solidVelocity, zero, zero, Zero);

				objectNodes.push_back(ip);
			}
		}
	}

	cout << objects.size() << " object(s) initialized in " << objectNodes.size() << " nodes\n";

}

void LB::updateObjectBoundaries(objectList& objects){



	findNewActiveObject(objects);

	const double zero = 0.0;
  tVect Zero = tVect(0.0,0.0,0.0);

	objectNodes.clear();

    for (int io = 0; io < objects.size(); io++) {
        const tVect convertedPosition = objects[io].x0 / unit.Length;
        const double convertedRadius = objects[io].r / unit.Length;
        const unsigned int indexHere = objects[io].index;
#pragma omp parallel for
        for (int it = 0; it < totNodes; ++it) {
            const tVect nodePosition = getPosition(it);
            if (nodePosition.insideSphere(convertedPosition, convertedRadius)) {// setting solidIndex
                types[it].setSolidIndex(indexHere);
                types[it].setObjectWall();
            }
        }
    }

	for (int ip = 0; ip < totNodes; ip++ ){
		if ( types[ip].isObjectWall() ) {

			delete nodes[ip];
			createNode(ip);
        	tVect solidVelocity;
        	const tVect nodePosition = getPosition(ip);
        	unsigned int solidIndex = types[ip].getSolidIndex();
		if (nodePosition.insideSphere(objects[solidIndex].x0 / unit.Length, objects[solidIndex].r / unit.Length)) {
            	solidVelocity = objects[solidIndex].x1 / unit.Speed;

        // reset velocity and mass (useful for plotting)
        // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
        		nodes[ip]->initialize(fluidMaterial.initDensity, solidVelocity, zero, zero, Zero);

				objectNodes.push_back(ip);
		}
		}
	}


}


void LB::findNewActiveObject(objectList& objects) {

    // SOLID TO ACTIVE CHECK

	//cout << "find new active object" << endl;
  tVect Zero = tVect(0.0,0.0,0.0);
    unsigned int maxX(0), maxY(0), maxZ(0);
    for (int it = 0; it < totNodes; ++it) {
        if (types[it].isActive()) {
            maxX = std::max(maxX, getX(it));
            maxY = std::max(maxY, getY(it));
            maxZ = std::max(maxZ, getZ(it));
        }
    }

    const tVect maxP(maxX, maxY, maxZ);

    for (int ip = 0; ip < objectNodes.size() ; ip++) {

//	if (types[ip].isObjectWall() && !types[ip].isInsideParticle){
	    const unsigned int index = objectNodes[ip];

        const tVect nodePosition = getPosition(index);
        // solid index to identify cluster
	   	const unsigned int objectIndex = 0; //types[index].getSolidIndex();

		tVect solidVelocity = objects[objectIndex].x1 / unit.Speed;
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        bool newActive = true;

	   	bool surroundedGas = true;
	       for (int j = 1; j < lbmDirec; ++j) {
	            const unsigned int link = nodes[index]->d[j];
	            if (types[link].isFluid()) {
	                surroundedGas = false;
	                break;
	            }
	        }

        if (nodePosition.insideSphere(objects[objectIndex].x0 / unit.Length, objects[objectIndex].r / unit.Length))
			{
                	newActive = false;
				//	if (!types[index].isInsideParticle()){
					//	delete nodes[index];
				//	}
			  }





		if (newActive && surroundedGas) {

			types[index].setGas();
			nodes[index] = 0;

			newActive = false;

		}




		if (newActive) {

			delete nodes[index];

			createNode(index);

			types[index].setFluid();

 	       for (int j = 1; j < lbmDirec; ++j) {
 	            const unsigned int link = nodes[index]->d[j];
 	            if (types[link].isGas()) {
 	                types[index].setInterface();
 	                break;
 	            }
 	        }

			if (types[index].isFluid()){
	            tVect deltah = getPosition(index) - maxP;
	            // density=initDensity; velocity=initVelocity, mass=initDensity; viscosity=initVisc; force=lbF
				nodes[index]->initialize(fluidMaterial.initDensity, initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF);
				//nodes[index]->initialize(fluidMaterial.initDensity + 3.0 * fluidMaterial.initDensity * deltah.dot(lbF), initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF);

			}
			else if (types[index].isInterface()){

			nodes[index]->initialize(fluidMaterial.initDensity, Zero, 0.5 * fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF);

			}



		}

	}
}


void LB::initializeTopography() {

    // SOLID WALLS ////////////////////////
    if (lbTopography) {
        cout << "Setting topography" << endl;
        lbTop.readFromFile(lbTopographyFile);
        lbTop.show();
        // check if topography grid contains the fluid domain
        ASSERT(lbTop.coordX[0] < 0.0);
        ASSERT(lbTop.coordY[0] < 0.0);
        ASSERT(lbTop.coordX[lbTop.sizeX - 1] > lbSize[0] * unit.Length);
        ASSERT(lbTop.coordY[lbTop.sizeY - 1] > lbSize[1] * unit.Length);
#pragma omp parallel for
        for (int it = 0; it < totNodes; ++it) {
            // control is done in real coordinates
            const tVect nodePosition = getPosition(it) * unit.Length;
            const double distanceFromTopography = lbTop.distance(nodePosition);
            //cout<<distanceFromTopography<<endl;
            if (distanceFromTopography < 0.0) {// setting solidIndex
                types[it].setTopography();
            }
        }


    }
}


void LB::initializeInterface() {
    // creates an interface electing interface cells from active cells
    // defining a surface

    cout << "Initializing gas nodes" << endl;
	unsigned int typeHere;
	double convertedFluidHeight = fluidHeight / unit.Length;
	for (int it = 0; it < totNodes; ++it) {
		if (getY(it) > convertedFluidHeight ) {
			typeHere = 2 ; // gas
			if (types[it].isFluid()) {   //not solid
				types[it].setType(typeHere);
			}
		}
	}

	cout << "Gas nodes initialized" << endl;

}

void LB::restartInterface(ifstream& fluidFileID) {
    // creates an interface electing interface cells from active cells
    // defining a surface

    cout << "Restarting interface" << endl;

    // checking for boundary between gas and fluid and assigning interface properties
    unsigned int typeHere;
    for (int it = 0; it < totNodes; ++it) {
        fluidFileID>>typeHere;
        if (types[it].isFluid()) {
            types[it].setType(typeHere);
        }
    }

}

void LB::initializeLists() {

	 //  cout << "Resetting lists" << endl;

    // note that interface is not defined here. All fluid, interface and gas cells are 0 at the moment
    fluidNodes.clear();
    particleNodes.clear();
    interfaceNodes.clear();

    // creating list and initialize macroscopic variables for all nodes except walls
    for (int it = 0; it < totNodes; ++it) {

        // FLUID NODES ////
        if (types[it].isFluid()) {
            // adding the free node indexes to the list
            fluidNodes.push_back(it);
        }            // INTERFACE NODES ////
        else if (types[it].isInterface()) {
            // interface cells list update
            interfaceNodes.push_back(it);
        }

        // PARTICLE NODES ////
        if (types[it].isInsideParticle()) {
            // adding the particle node to the list
            particleNodes.push_back(it);
        }
    }

}

void LB::initializeVariables() {

//    cout << "Initializing variables" << endl;
    // note that interface is not defined here. All fluid, interface and gas cells are uninitialized at the moment
    // calculate maximum height of the fluid

    unsigned int maxX(0), maxY(0), maxZ(0);
    for (int it = 0; it < totNodes; ++it) {
        if (types[it].isActive()) {
            maxX = std::max(maxX, getX(it));
            maxY = std::max(maxY, getY(it));
            maxZ = std::max(maxZ, getZ(it));
        }
    }
    const tVect maxP(maxX, maxY, maxZ);

    // checking for boundary between gas and fluid and assigning interface properties
    // at this point fluid cells contain actual fluid cells and potential interface cells, so we create the node anyway
    double massFluid = 0.0;
    double massInterface = 0.0;
    for (int it = 0; it < totNodes; ++it) {
        if (types[it].isFluid()) {
            // creating node
            createNode(it);
            // check if it is interface
            for (int j = 1; j < lbmDirec; ++j) {
                unsigned int link = nodes[it]->d[j];
                if (types[link].isGas()) {
                    types[it].setInterface();
                    break;
                }
            }
        }
        // now assign macroscopic quantities accordingly
        // FLUID NODES ////
        if (types[it].isFluid()) {
            massFluid += 1.0;
            // setting macroscopic variables
            // density is calculated using hydrostatic profile
            tVect deltah = getPosition(it) - maxP;
            // density=initDensity; velocity=initVelocity, mass=initDensity; viscosity=initVisc; force=lbF
            nodes[it]->initialize(fluidMaterial.initDensity, initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF);
  //          nodes[it]->initialize(fluidMaterial.initDensity + 3.0 * fluidMaterial.initDensity * deltah.dot(lbF), initVelocity, fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF); // density
        }// INTERFACE NODES ////
        else if (types[it].isInterface()) {
            massInterface += 0.5;
            // setting macroscopic variables
            // density=initDensity; velocity=initVelocity, mass=0.5*initDensity; viscosity=initVisc; force=lbF
            nodes[it]->initialize(fluidMaterial.initDensity, initVelocity, 0.5 * fluidMaterial.initDensity, fluidMaterial.initDynVisc, lbF);
        }
    }
    cout << "Approximate volume = " << massFluid*unit.Volume << " (fluid body), " << massInterface*unit.Volume << " (interface), " << (massFluid + massInterface)*unit.Volume << " (tot), " << endl;


}

void LB::restartVariables(ifstream& fluidFileID) {

    cout << "Restarting variables" << endl;

    // creating list and initialize macroscopic variables for all nodes except walls
    double nHere = 0.0;
    double uXHere = 0.0;
    double uYHere = 0.0;
    double uZHere = 0.0;
    double massHere = 0.0;
    double viscHere = 0.0;
    double fHere[lbmDirec];
    for (int j = 1; j < lbmDirec; ++j) {
        fHere[j] = 0.0;
    }

    for (int it = 0; it < totNodes; ++it) {
        if (types[it].isActive()) {
            fluidFileID>>nHere;
            fluidFileID>>uXHere;
            fluidFileID>>uYHere;
            fluidFileID>>uZHere;
            const tVect uHere(uXHere, uYHere, uZHere);
            fluidFileID>>massHere;
            fluidFileID>>viscHere;
            for (int j = 0; j < lbmDirec; ++j) {
                fluidFileID >> fHere[j];
            }
            // creating node
            createNode(it);
            // setting macroscopic variables
            if (types[it].isFluid()) {
                nodes[it]->restart(nHere + fluidMaterial.initDensity, uHere, massHere + fluidMaterial.initDensity, viscHere, fHere);
            }
            else if (types[it].isInterface()) {
                // density=initDensity; velocity=initVelocity, mass=0.5*initDensity; viscosity=initVisc; force=lbF
                nodes[it]->restart(nHere + fluidMaterial.initDensity, uHere, massHere + fluidMaterial.initDensity, viscHere, fHere);
            }
        }
    }
}

void LB::initializeWalls(wallList& walls) {
	cout << "Initializing wall nodes" << endl;

	const double zero = 0.0;
  tVect Zero = tVect(0.0,0.0,0.0);
    wallNodes.clear();

    // initializing wall nodes
    // note that, in the hypothesis that these walls are not evolving, only nodes at the interface need creation
    for (int it = 0; it < totNodes; ++it) {
	   if (types[it].isWall()) {

           // creating node
           createNode(it);
           // initialize node
           // STATIC WALL NODES ////
           if (types[it].isStatWall() || types[it].isSlipStatWall() || types[it].isTopography()) {
               // reset velocity and mass (useful for plotting)
               // density=0.0; velocity=(0.0,0.0,0.0), mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
               nodes[it]->initialize(fluidMaterial.initDensity, Zero, zero, zero, Zero);
           }                    // DYNAMIC WALL NODES ////

		   else if (types[it].isDynWall() || types[it].isSlipDynWall() || types[it].isCylinderWall() ) {
               // need to define velocity. It could be part of a cylinder or wall, we check both
               tVect solidVelocity;
               const tVect nodePosition = getPosition(it);
               unsigned int solidIndex = types[it].getSolidIndex();
               // objects
               // reset velocity and mass (useful for plotting)
               // density=0.0; velocity=solidVelocity, mass=0.0; viscosity=0.0; force=(0.0,0.0,0.0)
               nodes[it]->initialize(fluidMaterial.initDensity, solidVelocity, zero, zero, Zero);
           }
           // add node to list
           wallNodes.push_back(it);
        }
    }
}


// integration functions

void LB::cleanLists() {
    // sorts the active-node list and removes duplicates

    std::sort(fluidNodes.begin(), fluidNodes.end());
    std::sort(interfaceNodes.begin(), interfaceNodes.end());
    std::sort(particleNodes.begin(), particleNodes.end());
	std::sort(objectNodes.begin(), objectNodes.end());

    for (int ind = objectNodes.size() - 1; ind > 0; --ind) {
        const unsigned int i = objectNodes[ind];
        const unsigned int j = objectNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-object!" << endl;
            objectNodes.erase(objectNodes.begin() + ind);
        }
    }


    for (int ind = fluidNodes.size() - 1; ind > 0; --ind) {
        const unsigned int i = fluidNodes[ind];
        const unsigned int j = fluidNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-fluid!" << endl;
            fluidNodes.erase(fluidNodes.begin() + ind);
        }
    }
    for (int ind = interfaceNodes.size() - 1; ind > 0; --ind) {
        const unsigned int i = interfaceNodes[ind];
        const unsigned int j = interfaceNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-interface!" << endl;
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }
    }
    for (int ind = particleNodes.size() - 1; ind > 0; --ind) {
        const unsigned int i = particleNodes[ind];
        const unsigned int j = particleNodes[ind - 1];
        if (i == j) {
            cout << "duplicate-particleNode!" << endl;
            particleNodes.erase(particleNodes.begin() + ind);
        }
    }

    // list with active nodes i.e. nodes where collision and streaming are solved
    // solid nodes, particle nodes and gas nodes are excluded
    activeNodes.clear();
    activeNodes.reserve(fluidNodes.size() + interfaceNodes.size());
    activeNodes.insert(activeNodes.end(), fluidNodes.begin(), fluidNodes.end());
    activeNodes.insert(activeNodes.end(), interfaceNodes.begin(), interfaceNodes.end());

    std::sort(activeNodes.begin(), activeNodes.end());

    for (int ind = activeNodes.size() - 1; ind > 0; --ind) {
        unsigned int i = activeNodes[ind];
        unsigned int j = activeNodes[ind - 1];
        if (types[i].isGas()) {
            activeNodes.erase(activeNodes.begin() + ind);
            cout << "wrong-gas!" << endl;
        }
        else if (types[i].isObjectWall()) {
            activeNodes.erase(activeNodes.begin() + ind);
            cout << "wrong-object!" << endl;
        }
		 else if (nodes[i]==0) {
            activeNodes.erase(activeNodes.begin() + ind);
            cout << "memory error - uninitialized active node!" << endl;
        } else if (i == j) {
            activeNodes.erase(activeNodes.begin() + ind);
            cout << "duplicate-active!" << endl;
        }
    }

    for (int it = 0; it < activeNodes.size(); ++it) {
        if (types[activeNodes[it]].isGas()) {
            cout << "GAS ERROR" << endl;
            exit(0);
        }
    }
}

void LB::reconstruction(particleList& particles) {
    // reconstruction of macroscopic variables from microscopic distribution
    // this step is necessary to proceed to the collision step


#pragma omp parallel for
    for (int it = 0; it < activeNodes.size(); ++it) {
        nodes[activeNodes[it]]->reconstruct();
    }

	//#pragma omp parallel for
 /*   for (int ip = 0; ip < activeNodes.size(); ++ip) {
        unsigned int index = activeNodes[ip];


        if (types[index].isFluid() && types[index].isInsideParticle()) {

            const unsigned int particleIndex = types[index].getSolidIndex();

			nodes[index] -> yDisp = - particles[particleIndex].disp.y;
			nodes[index] -> rDisp = particles[particleIndex].disp_r;
			nodes[index] -> press = particles[particleIndex].pressure;
			nodes[index] -> n1 = particles[particleIndex].normal1;
			nodes[index] -> n2 = particles[particleIndex].normal2;
			nodes[index] -> yyStress = particles[particleIndex].STotal.m11;
		}

				if (types[index].isFluid() && !types[index].isInsideParticle()) {
			nodes[index] -> yDisp = 1.0e-6;
			nodes[index] -> rDisp = 1.0e-6;
			nodes[index] -> press = 1.0e-6;
			nodes[index] -> n1 = 1.0e-6 ;
			nodes[index] -> n2 = 1.0e-6 ;
			nodes[index] -> yyStress = 1.0e-6 ;

		}
	}	*/

}

void LB::collision() {

    if (!forceField) {
        lbF.reset();
    }

    //    const double amplitude=lbF.norm();

#pragma omp parallel for
    for (int it = 0; it < activeNodes.size(); ++it) {

        const unsigned int index = activeNodes[it];
        // equilibrium distributions
        double feq[lbmDirec];
        // shift velocity field to F/2
        nodes[index]->shiftVelocity(lbF);
        // compute equilibrium distributions
        nodes[index]->computeEquilibrium(feq);
        //        if (types[index].isInsideParticle()) {
        //            nodes[index]->resetViscosity(initDynVisc);
        //        }
        if (fluidMaterial.rheologyModel==BINGHAM || fluidMaterial.rheologyModel==FRICTIONAL || fluidMaterial.rheologyModel==VOELLMY || fluidMaterial.turbulenceOn) {
            // compute shear rate tensor, find invariant calculate viscosity (Bingham)
            nodes[index]->computeApparentViscosity(feq, fluidMaterial);
        }
        // compute new distributions
        nodes[index]->solveCollision(feq);
        // add force term to new distributions
        nodes[index]->addForce(lbF);
    }

}

void LB::streaming(objectList& objects, wallList& walls) {
    // STREAMING STEP
    // coefficients for equilibrium
    static double C1 = 3.0;
    static double C2 = 4.5;
    static double C3 = 1.5;
    // coefficient for free-surface
    static const double C2x2 = 9.0;
    static const double C3x2 = 3.0;
    // coefficient for slip conditions
    static const double S1 = slipCoefficient;
    static const double S2 = (1.0 - slipCoefficient);
    // creating list for collision function
    double staticPres[lbmDirec];
    for (int j = 0; j < lbmDirec; j++) {
        staticPres[j] = fluidMaterial.initDensity * coeff[j];
    }
    // coefficient for bounce-back
    static const double BBCoeff = 2.0 * 3.0;
    // extra mass due to bounce-back and moving walls
    double extraMass = 0.0;

    // initializing wall forces
    for (int iw = 0; iw < walls.size(); ++iw) {
        walls[iw].FHydro.reset();
    }
    // initializing object forces
    for (int io = 0; io < objects.size(); ++io) {
        objects[io].FHydro.reset();
		objects[io].MHydro.reset();
    }

    //  Saving in support variables f->fs
#pragma omp parallel for
    for (int it = 0; it < activeNodes.size(); ++it) {
        nodes[activeNodes[it]]->store();
    }

    //  Streaming
#pragma omp parallel for ordered reduction(+:extraMass)
    // cycling through active nodes
    for (int in = 0; in < activeNodes.size(); ++in) {

        // index of active node
        const unsigned int it = activeNodes[in];

        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // getting neighbor index
            const unsigned int link = nodes[it]->d[j];

            // if neighbor is normal fluid cell what follows is true
            if (types[link].isActive()) {
                // streaming is solved normally
                nodes[it]->f[opp[j]] = nodes[link]->fs[opp[j]];
            }
            else if (types[link].isGas()) {
                // additional variables for equilibrium f computation
                const double usq = nodes[it]->u.norm2();
                const double vuj = nodes[it]->u.dot(v[j]);
                //streaming with constant pressure interface
                nodes[it]->f[opp[j]] = -nodes[it]->fs[j] + coeff[j] * fluidMaterial.initDensity * (2.0 + C2x2 * (vuj * vuj) - C3x2 * usq);
            }
            else if (types[link].isTopography()) {
                // simple bounce-back
                nodes[it]->f[opp[j]] = nodes[it]->fs[j];
            }

//            }// for moving walls there is simple bounce-back with velocity correction
            else if (types[link].isDynWall()) {
                // getting the index of the wall to compute force in the right object
                const int solidIndex = types[link].getSolidIndex();
                // velocity of the wall
                const tVect vel = nodes[link]->u;
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes[it]->n * coeff[j] * vel.dot(v[j]); // mass!!!!!

                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes[it]->bounceBackForce(j, staticPres, BBi);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                {
                    walls[solidIndex].FHydro += BBforce;
                }
                nodes[it]->f[opp[j]] = nodes[it]->fs[j] - BBi;
                // adding the extra mass to the surplus
                extraMass += BBi * nodes[it]->mass;
            }                // for walls there is simple bounce-back
            else if (types[link].isStatWall()) {
                // getting the index of the wall to compute force in the right object
                const int solidIndex = types[link].getSolidIndex();
                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
                const tVect BBforce = nodes[it]->bounceBackForce(j, staticPres, 0.0);
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                {
                    walls[solidIndex].FHydro += BBforce;
                }
                nodes[it]->f[opp[j]] = nodes[it]->fs[j];
            }// for walls there is simple bounce-back


            else if (types[link].isObjectWall()) {
                // getting the index of the wall to compute force in the right object
                const int solidIndex = types[link].getSolidIndex();
                // static pressure is subtracted in order to correctly compute buoyancy for floating objects
			    const tVect vel = nodes[link]->u;
                // variation in Bounce-Back due to moving object
			    const double BBi = BBCoeff * nodes[it]->n * coeff[j] * vel.dot(v[j]);
                const tVect BBforce = nodes[it]->bounceBackForce(j, staticPres, BBi);
				const tVect rad = getPosition(link) - objects[solidIndex].x0 / unit.Length;
                // updating force and torque on the object (lattice units). This point is critical since many nodes update the force on the same object (lattice units)
#pragma omp ordered
                {
                    objects[solidIndex].FHydro += BBforce;
					objects[solidIndex].MHydro += rad.cross(BBforce);
                }
                nodes[it]->f[opp[j]] = nodes[it]->fs[j] - BBi;
				// adding the extra mass to the surplus
				extraMass += BBi * nodes[it]->mass;
            }	else if (types[link].isSlipStatWall()) {
                if (j > 6) {
                    bool active1 = false;
                    bool active2 = false;
                    const unsigned int nodeCheck1 = nodes[it]->d[slip1Check[j]];
                    const unsigned int nodeCheck2 = nodes[it]->d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1 = true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2 = true;
                    }
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes[it]->f[opp[j]] = S1 * nodes[nodeCheck1]->fs[slip1[j]] + S2 * nodes[it]->fs[j];

                    } else if (!active1 && active2) {
                        // second
                        nodes[it]->f[opp[j]] = S1 * nodes[nodeCheck2]->fs[slip2[j]] + S2 * nodes[it]->fs[j];
                    } else {
                        // standard BB
                        nodes[it]->f[opp[j]] = nodes[it]->fs[j];
                    }
                } else {
                    // standard BB
                    nodes[it]->f[opp[j]] = nodes[it]->fs[j];
                }
            } else if (types[link].isSlipDynWall()) {

                // velocity of the wall
                const tVect vel = nodes[link]->u;
                // variation in Bounce-Back due to moving object
                const double BBi = BBCoeff * nodes[it]->n * coeff[j] * vel.dot(v[j]);
                if (j > 6) {
                    bool active1 = false;
                    bool active2 = false;
                    const unsigned int nodeCheck1 = nodes[it]->d[slip1Check[j]];
                    const unsigned int nodeCheck2 = nodes[it]->d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1 = true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2 = true;
                    }
                    // given the environment, perform the right operation
                    if (active1 && !active2) {
                        // first
                        nodes[it]->f[opp[j]] = S1 * nodes[nodeCheck1]->fs[slip1[j]] + S2 * (nodes[it]->fs[j] - BBi);
                        // adding the extra mass to the surplus
                        extraMass += S2 * nodes[it]->mass*BBi;
                    } else if (!active1 && active2) {
                        // second
                        nodes[it]->f[opp[j]] = S1 * nodes[nodeCheck2]->fs[slip2[j]] + S2 * (nodes[it]->fs[j] - BBi);
                        // adding the extra mass to the surplus
                        extraMass += S2 * nodes[it]->mass*BBi;
                    } else {
                        // standard BB
                        nodes[it]->f[opp[j]] = nodes[it]->fs[j] - BBi;
                        // adding the extra mass to the surplus
                        extraMass += nodes[it]->mass*BBi;
                    }
                } else {
                    // standard BB
                    nodes[it]->f[opp[j]] = nodes[it]->fs[j] - BBi;
                    // adding the extra mass to the surplus
                    extraMass += nodes[it]->mass*BBi;
                }
            } else {
                cout << link << " " << getX(link) << " " << getY(link) << " " << getZ(link) << " " << types[link].getType() << " TYPE ERROR" << endl;
                exit(0);
            }
        }
    }


    // redistributing extra mass due to bounce back to interface cells
    redistributeMass(extraMass);


    for (int w = 0; w < walls.size(); ++w) {
        walls[w].FHydro *= unit.Force;
    }

    for (int o = 0; o < objects.size(); ++o) {
        objects[o].FHydro *= unit.Force;
		objects[o].MHydro *= unit.Torque;
    }
}

// free surface functions

void LB::updateMass() {
    // refer to the article of Miller or of Svec et al. for this

#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->newMass = nodes[interfaceNodes[it]]->mass;
    }

    // mass for interface nodes is regulated by the evolution equation
#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        const int index = interfaceNodes[it];
        // additional mass streaming to/from interface
        double deltaMass = 0.0;
        //cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // getting neighbor index
            const unsigned int link = nodes[index]->d[j];
            // average liquid fraction
            if (types[link].isInterface()) {
                // average liquid fraction
                const double averageMass = 0.5 * (nodes[link]->mass + nodes[index]->mass);
                deltaMass += averageMass * nodes[index]->massStream(j);
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            } else if (types[link].isFluid()) {
                const double averageMass = 1.0;
                deltaMass += averageMass * nodes[index]->massStream(j);
                //                dummyNode.f[opp[j]]=nodes[link].fs[opp[j]];
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            //} else if (types[link].isGas()) {
                // averageMass = 0.0;
                //                averageMass=1.0;//*nodes[index].mass;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
            } else if (types[link].isDynWall() || types[link].isObjectWall() ) {
                const double averageMass = 1.0 * nodes[index]->mass; //1.0*nodes[index].mass;//0.5*nodes[index].u.norm()/(nodes[link].u.norm()+nodes[index].u.norm());
                deltaMass += averageMass * nodes[index]->massStream(j);
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            } else if (types[link].isCylinderWall()) {
                const double averageMass = 1.0 * nodes[index]->mass; //1.0*nodes[index].mass;//0.5*nodes[index].u.norm()/(nodes[link].u.norm()+nodes[index].u.norm());
                deltaMass += averageMass * nodes[index]->massStream(j);
                //              double BBi=2.0*3.0*dummyNode.n*coeff[j]*vel.dot(v[j]);
                //              dummyNode.f[opp[j]]=dummyNode.fs[j]-BBi;
                //                deltaMass+=averageMass*(nodes[index].f[opp[j]]-nodes[index].fs[j]);
                //                cout<<"delta mass ="<<averageMass<<","<<deltaMass<<"\n";
            } else if (types[link].isSlipDynWall()) {
                if (j > 6) {
                    bool active1 = false;
                    bool active2 = false;
                    const unsigned int nodeCheck1 = nodes[index]->d[slip1Check[j]];
                    const unsigned int nodeCheck2 = nodes[index]->d[slip2Check[j]];
                    // check for the environment
                    if (types[nodeCheck1].isActive()) {
                        active1 = true;
                    }
                    if (types[nodeCheck2].isActive()) {
                        active2 = true;
                    }
                    // given the environment, perform the right operation
                    double averageMass=0.0;
                    if (active1 && !active2) {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * (1.0 - slipCoefficient) * nodes[index]->mass;
                    } else if (!active1 && active2) {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * (1.0 - slipCoefficient) * nodes[index]->mass;
                    } else {
                        // adding the extra mass to the surplus
                        averageMass += 1.0 * nodes[index]->mass;
                    }
                    deltaMass += averageMass * nodes[index]->massStream(j);
                } else {
                    // adding the extra mass to the surplus
                    const double averageMass = 1.0 * nodes[index]->mass;
                    deltaMass += averageMass * nodes[index]->massStream(j);
                }
            }
        }
        nodes[index]->newMass += deltaMass;
    }
    // mass for fluid nodes is equal to density
#pragma omp parallel for
    for (int it = 0; it < fluidNodes.size(); ++it) {
        nodes[fluidNodes[it]]->mass = nodes[fluidNodes[it]]->n;
    }
#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->mass = nodes[interfaceNodes[it]]->newMass;
    }
}

void LB::updateInterface() {
    // updates the location of the interface, creating and deleting interface nodes

    // variable for storage of mass surplus
    double massSurplus = 0.0;

    // lists for "mutant" nodes
    unsIntList filledNodes, emptiedNodes;
    // interface nodes created due to interface evolution
    unsIntList newInterfaceNodes;

    // filling lists of mutant nodes and changing their type
    findInterfaceMutants(filledNodes, emptiedNodes);

    // fixing the interface (always one interface between fluid and gas)
    smoothenInterface(filledNodes, emptiedNodes, newInterfaceNodes, massSurplus);

    // updating characteristics of mutant nodes
    updateMutants(filledNodes, emptiedNodes, newInterfaceNodes, massSurplus);

    // remove isolated interface cells (both surrounded by gas and by fluid)
    removeIsolated(massSurplus);

    // distributing surplus to interface cells
    redistributeMass(massSurplus);

}

void LB::findInterfaceMutants(unsIntList& filledNodes, unsIntList& emptiedNodes) {
    // checks for nodes going through the following transitions
    // interface -> fluid (filled node)
    // interface -> gas (emptied node)
    filledNodes.clear();
    emptiedNodes.clear();

    // CHECKING FOR MUTANT NODES (old backward style for cycle...)
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        const unsigned int i = interfaceNodes[ind];
        // CHECKING FOR NEW FLUID NODES from filling
        if (nodes[i]->mass > nodes[i]->n) {
            // updating mutant list
            filledNodes.push_back(i);
            // updating type
            types[i].setFluid();
            // updating lists (interface -> fluid)
            fluidNodes.push_back(i);
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }// CHECKING FOR NEW GAS NODES from emptying
        else if (nodes[i]->mass < 0.0) {
            // updating mutant list
            emptiedNodes.push_back(i);
            // updating type
            types[i].setGas();
            // updating lists (interface ->noList)
            interfaceNodes.erase(interfaceNodes.begin() + ind);
        }
    }
}

void LB::smoothenInterface(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus) {

    newInterfaceNodes.clear();

    // CHECKING FOR NEW INTERFACE NODES from neighboring a new fluid node
    for (int it = 0; it < filledNodes.size(); ++it) {
        const unsigned int index = filledNodes[it];
        // cycling through neighbors
        for (int j = 1; j < lbmDirec; ++j) {
            // index of neighbor node
            const unsigned int link = nodes[index]->d[j];
            // checking if node is gas (so to be tranformed into interface)
            if (types[link].isGas()) {
                types[link].setInterface();
                newInterfaceNodes.push_back(link);
                // creating node
                createNode(link);
                // node is becoming active and needs to be initialized
                // same density and velocity; 1% of the mass
                nodes[link]->initialize(fluidMaterial.initDensity, nodes[index]->u, 0.01 * fluidMaterial.initDensity, nodes[index]->visc, nodes[index]->hydroForce + lbF);
                // the 1% of the mass is taken form the surplus
                massSurplus -= 0.01 * fluidMaterial.initDensity;
                //                activeNodes.push_back(nodes[link].index);
            }
        }
    }

    // CHECKING FOR NEW INTERFACE NODES from neighboring a new gas node
    for (int it = 0; it < emptiedNodes.size(); ++it) {
        const unsigned int index = emptiedNodes[it];
        for (int j = 0; j < lbmDirec; ++j) {
            const unsigned int link = nodes[index]->d[j];
            if (types[link].isFluid()) {
                types[link].setInterface();
                newInterfaceNodes.push_back(link);
                // characteristics are inherited by previous fluid cell. Only mass must be updated to 99% of initial mass
                nodes[link]->mass = 0.99 * nodes[link]->n;
                // the remaining 1% of the mass is added to the surplus
                massSurplus += 0.01 * nodes[link]->n;
                // removing from fluid nodes
                for (int ind = fluidNodes.size() - 1; ind >= 0; --ind) {
                    const unsigned int i = fluidNodes[ind];
                    if (i == link) {
                        fluidNodes.erase(fluidNodes.begin() + ind);
                    }
                }
            }
        }
    }
}

void LB::updateMutants(unsIntList& filledNodes, unsIntList& emptiedNodes, unsIntList& newInterfaceNodes, double& massSurplus) {

    // resetting new gas macroscopic quantities
    for (int it = 0; it < emptiedNodes.size(); ++it) {
        bool isNotInterface = true;
        for (int it1 = 0; it1 < newInterfaceNodes.size(); ++it1) {
            if (emptiedNodes[it] == newInterfaceNodes[it1]) {
                isNotInterface = false;
            }
        }
        if (isNotInterface) {
            // updating mass surplus
            massSurplus += nodes[emptiedNodes[it]]->mass; // maybe problems
            // deleting node
            delete nodes[emptiedNodes[it]];
            nodes[emptiedNodes[it]] = 0;
        }
    }

    // resetting new fluid macroscopic quantities
    for (int it = 0; it < filledNodes.size(); ++it) {
        bool isNotInterface = true;
        for (int it1 = 0; it1 < newInterfaceNodes.size(); ++it1) {
            if (filledNodes[it] == newInterfaceNodes[it1]) {
                isNotInterface = false;
            }
        }
        if (isNotInterface) {
            // updating mass surplus
            massSurplus += nodes[filledNodes[it]]->mass - nodes[filledNodes[it]]->n;
            // setting liquid fraction for new fluid cell (other macroscopic characteristics stay the same)
            nodes[filledNodes[it]]->mass = nodes[filledNodes[it]]->n;
        }
    }
    // resetting neighbors for newly created interface cells
    for (int it = 0; it < newInterfaceNodes.size(); ++it) {
        interfaceNodes.push_back(newInterfaceNodes[it]);
    }

}

void LB::removeIsolated(double& massSurplus) {
    // remove isolated interface cells (surrounded either by only fluid or only solid cells)

    // checking if it is surrounded by fluid (in that case is converted to fluid). Solid is an exception
    // reverse cycle is needed because of deletion function
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        const unsigned int i = interfaceNodes[ind];
        bool surroundedFluid = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const unsigned int link = nodes[i]->d[j];
            if (types[link].isGas()) {
                surroundedFluid = false;
                break;
            }
        }
        if (surroundedFluid) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW FLUID NODE from surrounding\n";
            // update mass storage for balance
            massSurplus += nodes[i]->mass - nodes[i]->n;
            // update characteristics (inherited from the gas node)
            nodes[i]->mass = nodes[i]->n;
            types[i].setFluid();
            interfaceNodes.erase(interfaceNodes.begin() + ind);
            fluidNodes.push_back(i);
        }
    }

    // checking if it is surrounded by gas (in that case is converted to gas)
    // or, better, if it is not connected to fluid (could be connected to walls or particles)
    for (int ind = interfaceNodes.size() - 1; ind >= 0; --ind) {
        const unsigned int i = interfaceNodes[ind];
        bool surroundedGas = true;
        for (int j = 1; j < lbmDirec; ++j) {
            const unsigned int link = nodes[i]->d[j];
            if (types[link].isFluid()) {
                surroundedGas = false;
                break;
            }
        }
        // updating mass surplus
        if (surroundedGas) {
            //            cout<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].z<<" NEW GAS NODE from surrounding\n";
            // update mass
            massSurplus += nodes[i]->mass;
            interfaceNodes.erase(interfaceNodes.begin() + ind);
            types[i].setGas();
            delete nodes[i];
            nodes[i] = 0;
        }
    }
}

void LB::redistributeMass(const double& massSurplus) {
    // redistribute the mass surplus among interface cells
    const double addMass = massSurplus / interfaceNodes.size();

#pragma omp parallel for
    for (int it = 0; it < interfaceNodes.size(); ++it) {
        nodes[interfaceNodes[it]]->mass += addMass;
    }
}

void LB::enforceMassConservation() {
    // calculate total mass
    double thisMass = 0.0;
    for (int ia = 0; ia < activeNodes.size(); ++ia) {
        const unsigned int index=activeNodes[ia];
        if (!types[index].isInsideParticle()) {
            thisMass += nodes[index]->mass;
        }
    }

    // mass deficit
    const double massDeficit = (thisMass - totalMass);
    // cout<<endl<<"This="<<thisMass<<" tot="<<totalMass;

    // fix it
    redistributeMass(-0.01 * massDeficit);

}

// particle coupling functions

void LB::updateIndices(particleList& particles) {
    // checks if solid indices need to be changed due to periodic shifting
    // this is a compulsory step when working with ghost particles


    initializeParticleBoundaries(particles);

    initializeLists();


    ////    cout<<"Updating indices\n";
    ////    #pragma omp parallel for
    //    for (int ip=0; ip<particleNodes.size(); ++ip) {
    //        const unsigned int index=particleNodes[ip];
    //        const tVect nodePosition=positions[index];
    //        for (int p=0; p<particles.size();++p) {
    //            if (nodePosition.insideSphere(particles[p].x0/unit.Length,particles[p].r/unit.Length)) { //-0.5?
    //                // if the node is inside the element, then set the indices again and leave
    //                types[index].setSolidIndex(particles[p].particleIndex);
    //                break;
    //            }
    //        }
    //    }

}

void LB::computeHydroForces(particleList& particles) {

    // initializing the elements forces (lattice units)
#pragma omp parallel for
    for (int n = 0; n < particles.size(); ++n) {
        //initializing this time step hydrodynamic force
        particles[n].FHydro = tVect(0.0, 0.0, 0.0);
        particles[n].MHydro = tVect(0.0, 0.0, 0.0);
        // initializing the fluid mass for buoyancy
        particles[n].fluidVolume = 0.0;

		particles[n].SHydro.reset();
    }



#pragma omp parallel for
    // cycling through active nodes
    for (int ip = 0; ip < activeNodes.size(); ++ip) {
        unsigned int index = activeNodes[ip];
        // resetting hydrodynamic forces on nodes
        nodes[index]->hydroForce.reset();

        if (types[index].isInsideParticle()) { // && types[index].isFluid()
            // getting the index of the particle to compute force in the right object
            const unsigned int particleIndex = types[index].getSolidIndex();
          //  const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
            // calculating velocity of the solid boundary at the node (due to rotation of particles)
            // vectorized radius (real units)
            const tVect radius = getPosition(index) - particles[particleIndex].x0 / unit.Length;
            // update velocity of the particle node (u=v_center+omega x radius) (real units)
            const tVect localVel = particles[particleIndex].x1 / unit.Speed + (particles[particleIndex].w1.cross(radius)) / unit.AngVel;

            // calculate differential velocity
            const tVect diffVel = nodes[index]->liquidFraction()*(nodes[index]->u - localVel);



            // force on fluid
            nodes[index]->hydroForce += -1.0 * diffVel; // -1.0*nodes[index]->liquidFraction()*diffVel;


            // force on particle
#pragma omp critical
            {

				particles[particleIndex].fluidVolume += nodes[index]->mass;
                particles[particleIndex].FHydro += 1.0 * diffVel;
                particles[particleIndex].MHydro += 1.0 * radius.cross(diffVel);

				particles[particleIndex].SHydro.m00 += 1.0 * radius.x * unit.Length * diffVel.x * unit.Force;
				particles[particleIndex].SHydro.m01 += 1.0 * radius.x * unit.Length * diffVel.y * unit.Force;
				particles[particleIndex].SHydro.m02 += 1.0 * radius.x * unit.Length * diffVel.z * unit.Force;
				particles[particleIndex].SHydro.m10 += 1.0 * radius.y * unit.Length * diffVel.x * unit.Force;
				particles[particleIndex].SHydro.m11 += 1.0 * radius.y * unit.Length * diffVel.y * unit.Force;
				particles[particleIndex].SHydro.m12 += 1.0 * radius.y * unit.Length * diffVel.z * unit.Force;
				particles[particleIndex].SHydro.m20 += 1.0 * radius.z * unit.Length * diffVel.x * unit.Force;
				particles[particleIndex].SHydro.m21 += 1.0 * radius.z * unit.Length * diffVel.y * unit.Force;
				particles[particleIndex].SHydro.m22 += 1.0 * radius.z * unit.Length *  diffVel.z * unit.Force;
            }
        }

    }


    // shifting elements forces and torques to physical units
#pragma omp parallel for
    for (int n = 0; n < particles.size(); ++n) {

        // adding buoyancy
        const tVect buoyancy = (-particles[n].dens) * particles[n].fluidVolume*lbF;
	  //    particles[n].FHydro += buoyancy;
        particles[n].FHydro *= unit.Force;
        particles[n].MHydro *= unit.Torque;
        particles[n].fluidVolume *= unit.Volume;

    }
}

void LB::findNewActive(unsIntList& newPopUpNodes, particleList& particles) {

    // SOLID TO ACTIVE CHECK
    // cycling through particle nodes
    //    #pragma omp parallel for ordered
    for (int ip = 0; ip < particleNodes.size(); ++ip) {
        const unsigned int index = particleNodes[ip];
        const tVect nodePosition = getPosition(index);
        // solid index to identify cluster
        const unsigned int particleIndex = types[index].getSolidIndex();
   //     const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
        // in this case check if it has been uncovered (must be out of all particles of the cluster) - we start with a true hypothesis
        bool newActive = true;
        // cycling through component particles
       // unsigned int componentIndex = 0;
      //  for (int j = 0; j < elmts[clusterIndex].components.size(); ++j) {
            // getting indexes from particle composing the cluster
         //   componentIndex = elmts[clusterIndex].components[j];
            // checking if it has been uncovered in component j of the cluster
            // radius need to be increased by half a lattice unit
            // this is because solid boundaries are located halfway between soli and fluid nodes
        if (nodePosition.insideSphere(particles[particleIndex].x0 / unit.Length, particles[particleIndex].r / unit.Length))
			{
                // if the node is still inside the element, the hypothesis of new active is not true anymore
                	newActive = false;
                // and we can get out of the cycle
                break;
            }
       // }

		if (newActive) {
            // turning up the cell
            //            #pragma omp ordered
            newPopUpNodes.push_back(index);
            types[index].setOutsideParticle();
            //            cout<<"new active\n";
        }
    }

    for (int it = 0; it < newPopUpNodes.size(); ++it) {
        //        cout<<"New active\n";
        for (int ind = particleNodes.size() - 1; ind >= 0; --ind) {
            if (newPopUpNodes[it] == particleNodes[ind])
			{
                // deleting the cell from particle list
                particleNodes.erase(particleNodes.begin() + ind);
            }
        }
    }
}


void LB::findNewSolid(unsIntList& newSolidNodes, particleList& particles) {

    // ACTIVE TO SOLID CHECK
    // we check only first order neighbors of particle nodes. This is done in order to avoid cycling through all active cells
    //    #pragma omp parallel for ordered
    for (int it = 0; it < particleNodes.size(); ++it) {

        const unsigned int index = particleNodes[it];
        //    for (intlist::iterator it=particleNodes.begin(); it<particleNodes.end(); ++it) {
        // solid index to identify cluster
        const unsigned int particleIndex = types[index].getSolidIndex();
       // const unsigned int clusterIndex = particles[particleIndex].clusterIndex;
        // cycle through first neighbors
        for (int k = 1; k < lbmMainDirec; ++k) {
            const unsigned int link = nodes[index]->d[k];
            // checking if solid particle is close to an active one -> we have an active node to check
            if (!types[link].isInsideParticle()) {
                const tVect linkPosition = getPosition(link);
                // check if neighbors has been covered (by any of the particles of the cluster) - we start with a false hypothesis
                bool newSolid = false;
         //       unsigned int componentIndex = 0;
                // cycling through all components of the cluster
           //     for (int j = 0; j < elmts[clusterIndex].components.size(); ++j) {
                    //                    cout<<j<<"\n";
                    // getting component particle index
             //       componentIndex = elmts[clusterIndex].components[j];

				    // check if it getting inside
                    // radius need to be increased by half a lattice unit
                    // this is because solid boundaries are located halfway between soli and fluid nodes
          if (linkPosition.insideSphere(particles[particleIndex].x0 / unit.Length, particles[particleIndex].r / unit.Length))
		   		   { //-0.5?
                        // if so, then the false hypothesis does not hold true anymore
                       newSolid = true;
                        // and we exit the cycle
                       break;
                   }
             // }
                // an active cell can become part of multiple particle at the same time
                // we need to check if it is already on the list
                if (newSolid) {
                    // check if we are creating a duplicate - we start with a false hypothesis
                    bool alreadyInside = false;
                    for (int it1 = 0; it1 < newSolidNodes.size(); ++it1) {
                        //                    for (intlist::iterator it1=newSolidNodes.begin(); it1!=newSolidNodes.end(); ++it1) {
                        if (link == newSolidNodes[it1]) {
                            // if it is already in the list than the true hypothesis does not hold true anymore
                            alreadyInside = true;
                            // and we exit the cycle
                            break;
                        }
                    }
                    // at this point, if it is both newSolid and not alreadyInside we can add it to particle nodes
                    if (!alreadyInside) {
                        // solid index is assigned here to avoid sending this information to the switching functions
                        //                        # pragma omp ordered
                        {
                            types[link].setSolidIndex(particleIndex);
                            types[link].setInsideParticle();
                            // sending node index for updating
                            newSolidNodes.push_back(link);
                            particleNodes.push_back(link);
                        }
                    }
                }
            }
        }
    }
}

void LB::solidToActive(unsIntList& newPopUpNodes,  double& massSurplus) {

}

void LB::activeToSolid(unsIntList& newSolidNodes, double& massSurplus) {


}

// functions for index management

tVect LB::getPosition(const unsigned int& index) const {
    unsigned int x, y, z;

    // index is calculated in this fashion:
    // index = x + y*X + z*X*Y
    // where X and Y are sizes of the lattice in x and y direction

    // from this stems that
    // x + y*X = index MOD X*Y
    // x = x + y*X MOD X
    // y = x + y*X DIV X
    // z = index DIV X*Y

    // see online documentation for class div_t (stdlib.h)
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    secondDiv = div(firstDiv.rem, int(lbSize[0]));

    x = secondDiv.rem;
    y = secondDiv.quot;
    z = firstDiv.quot;

    return tVect(double(x), double(y), double(z));
}

unsigned int LB::getX(const unsigned int& index) const {

    // see function getPosition for documentation
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    secondDiv = div(firstDiv.rem, int(lbSize[0]));

    return secondDiv.rem;
}

unsigned int LB::getY(const unsigned int& index) const {

    // see function getPosition for documentation
    div_t firstDiv, secondDiv;

    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));
    secondDiv = div(firstDiv.rem, int(lbSize[0]));

    return secondDiv.quot;
}

unsigned int LB::getZ(const unsigned int& index) const {

    // see function getPosition for documentation
    div_t firstDiv;

    firstDiv = div(int(index), int(lbSize[0] * lbSize[1]));

    return firstDiv.quot;
}
