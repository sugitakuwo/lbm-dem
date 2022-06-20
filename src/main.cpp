#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <csignal>
#include <limits>

#include "getpot.h"

#include "IO.h"
#include "DEM.h"
#include "LB.h"

/*

Code to simulate free-falling impactor onto dense suspensions

I wrote this code by modifying the Hybird DEM-LBM code by Alessandro Leonardi, changed the particle interactions to 3D-DEM by Satoshi Takada, and used the neigbor list and lubrication calculation of Susp3d by Anthony Ladd.

 */

enum ExitCode {
    UNFINISHED = -1, SUCCESS = 0, TIME_LIMIT_REACHED, SIGNAL_CAUGHT, ERROR
};

ExitCode exit_code = UNFINISHED;
ProblemName problemName = NONE;

void catchSignal(int sig) {
    std::cout << "Received signal " << sig << ": " << strsignal(sig) << " " << std::endl;
    exit_code = (SIGNAL_CAUGHT > exit_code ? SIGNAL_CAUGHT : exit_code);
}

void print_help() {
    cout << "Usage: " << " [-c configFile] [-d datadir] [-n resultFolderName (opt. time)] [options}\n";
}

void goCycle(IO& io, DEM& dem, LB& lb) {

     // advance one step in time
    io.realTime+=lb.unit.Time;
    ++io.currentTimeStep;
    ++lb.time;

if (io.demSolve) {
	    dem.discreteElementStep(lb.boundary);
	  }

if (io.lbmSolve && dem.demTime >= dem.demInitialRepeat) {
    if (io.demSolve) {
        lb.latticeBoltzmannCouplingStep(dem.newNeighborList, dem.particles, dem.objects);
    }

    lb.latticeBolzmannStep(dem.particles, dem.walls, dem.objects);

    if (lb.freeSurface) {
        lb.latticeBoltzmannFreeSurfaceStep();
    }
} else if (io.lbmSolve && dem.demTime <= dem.demInitialRepeat) {

    if (io.demSolve) {
        lb.latticeBoltzmannCouplingStep(dem.newNeighborList,dem.particles,dem.objects);
    }
}
   io.outputStep(lb, dem);

}

void parseCommandLine(IO& io, GetPot& commandLine) {

    //io.base_directory = command_line.follow("./", "--directory");

    // print help and exit if requested
    if (commandLine.search("-h") || commandLine.search("--help")) {
        print_help();
        exit(0);
    }

    // create the data directory in case it doesn't exist
    io.workDirectory = "data";
    if (commandLine.search("-d")) {
        io.workDirectory = commandLine.next(io.workDirectory);
    }

    // if the name of the simulation is specified create subfolder
    std::string programName = "name";
    if (commandLine.search("-n")) {
        // creating work directory
        programName = commandLine.next(programName);
        // if the specified name is "time" then create a folder named with the initial time
        if (programName == "time") {
            ostringstream convertTime;
            convertTime << std::setfill('0') << std::setw(4) << (io.now->tm_year + 1900)
                    << std::setfill('0') << std::setw(2) << (io.now->tm_mon + 1)
                    << std::setfill('0') << std::setw(2) << io.now->tm_mday ;
            string simulationTime = convertTime.str();
      //      io.workDirectory = io.workDirectory + "/" + simulationTime;
            io.workDirectory = io.workDirectory.substr(0, io.workDirectory.size());
        } else {
            io.workDirectory = io.workDirectory + "/" + programName;
        }

        //        io.workDirectory = command_line.next(io.workDirectory);
    }

    // lbm configuration file
    io.configFileName = "config.cfg";
    if (commandLine.search("-c")) {
        io.configFileName = commandLine.next(io.configFileName.c_str());
        cout << "Using " << io.configFileName << " for LBM and DEM parameters\n";
        // copy the lbm configuration file into the data folder for later reproducibility
        int a;
        a = mkdir(io.workDirectory.c_str(), 0777);
        cout << "Work directory created = " << io.workDirectory << ". Result: " << a << "\n";

        std::system(("cp '" + io.configFileName + "' '" + io.workDirectory + "'").c_str());
    }

    // make sure the config files can be read
    std::ifstream configFile(io.configFileName.c_str());
    if (!configFile) {
        cout << "ERROR: Can't open config file \"" << io.configFileName << "\" for reading!\n";
        //        return ERROR;
    }
    configFile.close();
}

void parseConfigFile(IO& io, DEM& dem, LB& lb, GetPot& configFile, GetPot& commandLine) {

    // GETTING SIMULATION PARAMETERS  /////////
    // DEM initial iterations
    PARSE_CLASS_MEMBER(configFile, dem.demInitialRepeat, "demInitialRepeat", 0.0);
    ASSERT(dem.demInitialRepeat >= 0);
    // LB iteration without coupling (for initial stability) - > time is frozen here
    PARSE_CLASS_MEMBER(configFile, lb.lbmInitialRepeat, "lbmInitialRepeat", 0);
    ASSERT(lb.lbmInitialRepeat >= 0);
    // maximum time variable value
    PARSE_CLASS_MEMBER(configFile, io.maxTime, "maxTime", 0.0);
    ASSERT(io.maxTime >= 0);

	cout << "maxTime= " << io.maxTime << endl;


    // for time integration and output
    PARSE_CLASS_MEMBER(configFile, io.maximumTimeSteps, "maximumTimeSteps", 0);
    ASSERT(io.maximumTimeSteps >= 0);
    PARSE_CLASS_MEMBER(configFile, io.screenExpTime, "screenExpTime", 0.0);
    ASSERT(io.screenExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidExpTime, "fluidExpTime", 0.0);
    ASSERT(io.fluidExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidLagrangianExpTime, "fluidLagrangianExpTime", 0.0);
    ASSERT(io.fluidLagrangianExpTime >= 0);
  //  PARSE_CLASS_MEMBER(configFile, io.fluid2DExpTime, "fluid2DExpTime", 0.0);
  //  ASSERT(io.fluid2DExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.partExpTime, "partExpTime", 0.0);
    ASSERT(io.partExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.fluidRecycleExpTime, "fluidRecycleExpTime", 0.0);
    ASSERT(io.fluidRecycleExpTime >= 0);
    PARSE_CLASS_MEMBER(configFile, io.partRecycleExpTime, "partRecycleExpTime", 0.0);
    ASSERT(io.partRecycleExpTime >= 0);

	io.velExpTime = io.chainExpTime = io.objectExpTime = io.objectExpTime2 = io.bodyExpTime = io.springExpTime = io.boundExpTime = io.partExpTime;

	io.simComplexExpTime = io.part_pressureExpTime = io.part_yyStressExpTime = io.part_n2ExpTime = io.part_n1ExpTime =io.part_rDispExpTime =io.part_yDispExpTime = io.partExpTime;

	io.stressExpTime = io.dispExpTime = io.fluidExpTime;


    // solving particles?
    PARSE_CLASS_MEMBER(configFile, io.demSolve, "demSolve", 0);
    // solving fluid?
    PARSE_CLASS_MEMBER(configFile, io.lbmSolve, "lbSolve", 0);
    // solving rotational reference system (= adding extra apparent accelerations?)

    PARSE_CLASS_MEMBER(configFile, lb.freeSurface, "freeSurfaceSolve", 0);
    // is there a force field? (shutting down improves performances)
    PARSE_CLASS_MEMBER(configFile, lb.forceField, "forceFieldSolve", 0);

    // getting LBM parameters
    lb.latticeBoltzmannGet(configFile, commandLine);
    // getting DEM parameters and particle initial status
    dem.discreteElementGet(configFile, commandLine);


}

void printUfo(GetPot& command_line, GetPot& configFile) {
    // warn about unused parameters
    std::vector<std::string> ufos = configFile.unidentified_variables();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized lbm config file parameter(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
    ufos = command_line.unidentified_arguments();
    if (ufos.size() > 0) {
        cout << INFO << "WARNING: Unrecognized command line argument(s):";
        for (unsigned int i = 0; i < ufos.size(); ++i) {
            cout << " " << ufos[i];
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {

    // Checking number of processes involved
    cout << "Program starts with threads:";
#pragma omp parallel
    {
#pragma omp critical
        cout << " X ";
    }
    cout << "\n";

    // DECLARATION OF VARIABLES - Input-Output ///////////////
    IO io;

    // DECLARATION OF VARIABLES - DEM ///////////////
    DEM dem;

    // DECLARATION OF VARIABLES - LB ///////////////
    LB lb;

    // print some info for restorability
    time_t t = time(0); // get time now
    io.now = localtime(&t);
    cout << "Binary was compiled on " << __DATE__ << " at " << __TIME__ << endl;
    cout << "Invocation on " << (io.now->tm_year + 1900) << "-" << (io.now->tm_mon + 1) << "-" << io.now->tm_mday <<
            ", at " << io.now->tm_hour << ":" << io.now->tm_min << ":" << io.now->tm_sec << endl;
    cout << "Invocation line: " << argv[0];
    for (int i = 1; i < argc; ++i) {
        cout << " " << argv[i];
    }
    cout << endl;

    // parsing command line and configuration file
    GetPot commandLine(argc, argv);
    parseCommandLine(io, commandLine);

    // parsing LBM input file
    GetPot configFile(io.configFileName);
    parseConfigFile(io, dem, lb, configFile, commandLine);

    printUfo(commandLine, configFile);

    // bind signal handlers
    signal(SIGHUP, SIG_IGN); // for easy running through ssh
    signal(SIGQUIT, catchSignal);
    signal(SIGTERM, catchSignal); // important for job killing on a cluster

    io.currentTimeStep = 0;

    /* /////////////////////////////////////////////////////////////////////////
    // PROGRAM CORE  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////*/

    //  GLOBAL INITIALIZATION ///////////////////////////////////////////////
    // initializing input-output
    io.initialize();

    // initializing lattice
    lb.latticeDefinition();
    lb.LBShow();

    // initializing DEM parameters
    const tVect externalForce = lb.lbF * lb.unit.Accel;
    const tVect externalRotation = lb.lbRot * lb.unit.AngVel;
    const tVect externalRotationCenter = lb.lbRotCenter * lb.unit.Length;

    dem.discreteElementInit(lb.boundary, lb.lbPhysicalSize, lb.lbBoundaryLocation, lb.unit.Time, externalForce, lb.fluidHeight, lb.freeSurface);

    if (io.lbmSolve) {
        lb.latticeBolzmannInit(dem.walls, dem.particles, dem.objects);
    }

    // setting time
    lb.time = 0;
    dem.demTime = 0.0;
    io.realTime=0.0;
    dem.demTimeStep = 0;

    // initial output
    io.outputStep(lb, dem);

    // CYCLE /////////////////////////////
    // integrate in time
    while (true) {

        if (io.realTime != 0.0 && io.realTime > io.maxTime) {
            exit_code = SUCCESS;
        }            // exit normally if the maximum simulation time has been reached
      		else {
            // core of the code, performs time steps
            goCycle(io, dem, lb);

            //            // exit abnormally if a serious problem has occurred
            //            if (io.problem) {
            //                exit_code = ERROR;
            //            }
        }

	  if (exit_code > UNFINISHED) {
		       break;
		  }
    }
    return exit_code;
}
