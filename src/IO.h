
#ifndef IO_H
#define IO_H


#include <cstdlib>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <ctime>
#include <csignal>
#include <iomanip>
// only to create a directory
#include <sys/stat.h>
#include <sys/types.h>
//
#include "LB.h"
#include "DEM.h"
#include "myvector.h"
#include "elmt.h"
#include "node.h"

using namespace std;
extern ProblemName problemName;

// define precision outputs for statistic file
#define wSt int(15)

class IO {

public:
    // simulation time in real units
    double realTime;
    // simulation iterations so far
    unsigned int currentTimeStep;
    // switchers for solvers
    bool demSolve, lbmSolve;
    // scwitchers for rotation of the reference frame
    // locations
    string workDirectory, partDirectory, fluidDirectory;
    // time at invocation
    struct tm *now;
    // config file
    string configFileName;
    // time integration: max number of iterations
    unsigned int maximumTimeSteps;
    // time integration: maximum time
    double maxTime;
    // intervals for output
    double screenExpTime, stateExpTime, fluidExpTime, fluidLagrangianExpTime, fluid2DExpTime, partExpTime, boundExpTime, objectExpTime, bodyExpTime, springExpTime, objectExpTime2, chainExpTime, velExpTime,dispExpTime,stressExpTime;
    double part_yDispExpTime, part_rDispExpTime, part_n1ExpTime, part_n2ExpTime, part_pressureExpTime, part_yyStressExpTime;
	// recycle intervals
    double fluidRecycleExpTime, partRecycleExpTime, outputExpTime;
	  double simComplexExpTime;

    // general functions
    void initialize();
    void outputStep(const LB& lb, DEM& dem); // DEM used to be passed as const
    void outputFinal();

private:

    // generic file output streams
    string exportFileName;
    ofstream exportFile;
    unsigned int lastScreenExp;


    // function that groups file creations
    void createFiles(const LB& lb, const DEM& dem);

    // RECYCLE /////////////////////////////////////////////////////////////////////////////////////////

    // particle recycle file
    unsigned int lastPartRecycleExp;
    string partRecycleFileFormat;
    void exportRecycleParticles(const particleList& particles, const string& partRecycleFile);

    // simplicial complex file
    unsigned int lastSimComplexExp;
    string simComplexFileFormat;
    void exportSimComplex(const DEM& dem, const string& simComplexFile);

    // fluid recycle file
    unsigned int lastFluidRecycleExp;
    string fluidRecycleFileFormat;
    void exportRecycleFluid(const LB& lb, const string& fluidRecycleFile);

    // PARAVIEW /////////////////////////////////////////////////////////////////////////////////////////

    // particle paraview file
    unsigned int lastPartExp;
    string partFileFormat;
    void exportParaviewParticles(const particleList& particles, const string& particleFile);

    unsigned int lastPart_yDispExp;
    string part_yDispFileFormat;
    void exportParaviewParticles_yDisp(const particleList& particles, const string& particle_yDispFile);

    unsigned int lastPart_rDispExp;
    string part_rDispFileFormat;
    void exportParaviewParticles_rDisp(const particleList& particles, const string& particle_rDispFile);

    unsigned int lastPart_n1Exp;
    string part_n1FileFormat;
    void exportParaviewParticles_n1(const particleList& particles, const string& particle_n1File);

    unsigned int lastPart_n2Exp;
    string part_n2FileFormat;
    void exportParaviewParticles_n2(const particleList& particles, const string& particle_n2File);

    unsigned int lastPart_pressureExp;
    string part_pressureFileFormat;
    void exportParaviewParticles_pressure(const particleList& particles, const string& particle_pressureFile);

    unsigned int lastPart_yyStressExp;
    string part_yyStressFileFormat;
    void exportParaviewParticles_yyStress(const particleList& particles, const string& particle_yyStressFile);

	// object paraview file
    unsigned int lastObjectExp;
    string objectFileFormat;
    void exportParaviewObjects(const objectList& objects, const string& objectFile);

    unsigned int lastBodyExp;
    string bodyFileFormat;
    void exportParaviewBodies(const bodyList& bodies, const string& bodyFile);

    unsigned int lastSpringExp;
    string springFileFormat;
    void exportParaviewSprings(const bodyList& bodies, const objectList& objects, const string& springFile);

    unsigned int lastObjectExp2;
    string objectFileFormat2;
    void exportParaviewObjects2(const objectList& objects, const string& objectFile2);

	// boundary paraview file
	unsigned int lastBoundExp;
    string boundFileFormat;
    void exportParaviewBoundaries(const DEM& dem, const string& boundaryFile);

	unsigned int lastVelExp;
    string velFileFormat;
    void exportParaviewParticleVelocities(const particleList& particles, const string& velFile);

	unsigned int lastChainExp;
    string forceChainsFileFormat;
    void exportForceChains(const DEM& dem, const string& forceChainsFile);




    // Eulerian fluid paraview file
    unsigned int lastFluidExp;
    string fluidFileFormat;
    void exportEulerianParaviewFluid(const LB& lb, const string& fluidFile);

    unsigned int lastDispExp;
    string dispFileFormat;
    void exportParaviewDisp(const LB& lb, const string& dispFile);

    unsigned int lastStressExp;
    string stressFileFormat;
    void exportParaviewStress(const LB& lb, const string& stressFile);

    // Lagrangian fluid paraview file
    unsigned int lastFluidLagrangianExp;
    string fluidLagrangianFileFormat;
    void exportLagrangianParaviewFluid(const LB& lb, const string& fluidFile);



    // GENERAL OUTPUT FILES ////////////////////////////////////////////////////////////////////////////////////

	string sedimentRateFileName;
	ofstream sedimentRateFile;
	void exportSedimentRate(const DEM& dem);

	string objectPositionFileName;
	ofstream objectPositionFile;
	void exportObjectPosition(const objectList& objects, const LB& lb);

  string bodyPositionFileName;
  ofstream bodyPositionFile;
  void exportBodyPosition(const bodyList& bodies);

	string objectForceFileName;
	ofstream objectForceFile;
	void exportObjectForce(const objectList& objects, const LB& lb);

  string objectPositionFileName2;
  ofstream objectPositionFile2;
  void exportObjectPosition2(const objectList& objects, const LB& lb);

  string objectForceFileName2;
  ofstream objectForceFile2;
  void exportObjectForce2(const objectList& objects, const LB& lb);

	string wallForceFileName;
	ofstream wallForceFile;
	void exportWallForce(const wallList& walls);

	string stressTensorFileName;
	ofstream stressTensorFile;
	void exportStressTensor(const DEM& dem);

	string normalStressContFileName;
	ofstream normalStressContFile;
	void exportNormalStressContributions(const DEM& dem);

	string shearStressContFileName;
	ofstream shearStressContFile;
	void exportShearStressContributions(const DEM& dem);



};

#endif /* IO_H */
