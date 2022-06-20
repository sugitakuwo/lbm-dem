
#include "IO.h"

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PUBLIC FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/


void IO::initialize() {

    partDirectory = workDirectory + "/particleData";
    fluidDirectory = workDirectory + "/fluidData";

    int a;
    if (demSolve) {
        a = mkdir(partDirectory.c_str(), 0777);
        cout << "Work directory created = " << partDirectory << ". Result: " << a << "\n";
    }
    if (lbmSolve) {
        a = mkdir(fluidDirectory.c_str(), 0777);
        cout << "Work directory created = " << fluidDirectory << ". Result: " << a << "\n";
    }

    // build data file names
    fluidFileFormat = fluidDirectory + "/fluid%010u.vti";
	  stressFileFormat = fluidDirectory + "/stress%010u.vti";
	  dispFileFormat = fluidDirectory + "/disp%010u.vti";
    fluidLagrangianFileFormat = fluidDirectory + "/fluidLag%010u.vtu";
    partFileFormat = partDirectory + "/part%04u.vtk";
	   forceChainsFileFormat = partDirectory + "/forceChains%04u.vtk";
	    boundFileFormat = partDirectory + "/bound%04u.vtk";
	     velFileFormat = partDirectory + "/vel%04u.vtk";
    fluidRecycleFileFormat = fluidDirectory + "/fluidRecycle%010u.dat";
    partRecycleFileFormat = partDirectory + "/partRecycle%04u.dat";
	   simComplexFileFormat = partDirectory + "/simComplex%03u.txt";
    objectFileFormat = partDirectory + "/object%04u.vtk";
    objectFileFormat2 = partDirectory + "/2_object%04u.vtk";
    bodyFileFormat = partDirectory + "/body%04u.vtk";
    springFileFormat = partDirectory + "/spring%04u.vtk";




    //  initializing output file
    exportFileName = workDirectory + "/export.dat";

	if (demSolve){


		objectPositionFileName = workDirectory + "/objPosition.dat";
    	objectPositionFile.open(objectPositionFileName.c_str(), ios::app);
    	objectPositionFile << "position(x) position(y) position(z) velocity(x) velocity(y) velocity(z) depth(y)  \n";
    	objectPositionFile.close();

        bodyPositionFileName = workDirectory + "/bodyPosition.dat";
        bodyPositionFile.open(bodyPositionFileName.c_str(), ios::app);
        bodyPositionFile << "position(x) position(y) position(z) velocity(x) velocity(y) velocity(z) Length Delta  \n";
        bodyPositionFile.close();

		objectForceFileName = workDirectory + "/objForce.dat";
    	objectForceFile.open(objectForceFileName.c_str(), ios::app);
    	objectForceFile << "FExerted(x) FExerted(y) FExerted(z) SExerted(x) SExerted(y) SExerted(z) FTotal(x) FTotal(y) FTotal(z) FContact(x) FContact(y) FContact(z) FHydro(x) FHydro(y) FHydro(z) FLub(x) FLub(y) FLub(z) FRep(x) FRep(y) FRep(z) FGrav(x) FGrav(y) FGrav(z) \n";
    	objectForceFile.close();



		wallForceFileName = workDirectory + "/wallForce.dat";
    	wallForceFile.open(wallForceFileName.c_str(), ios::app);
    	wallForceFile << "0-FTotal(x) 0-FTotal(y) 0-FTotal(z) 1-FTotal(x) 1-FTotal(y) 1-FTotal(z) 2-FTotal(x) 2-FTotal(y) 2-FTotal(z) 3-FTotal(x) 3-FTotal(y) 3-FTotal(z) 4-FTotal(x) 4-FTotal(y) 4-FTotal(z) 5-FTotal(x) 5-FTotal(y) 5-FTotal(z) \n";
    	wallForceFile.close();

		stressTensorFileName = workDirectory + "/stressTensor.dat";
    	stressTensorFile.open(stressTensorFileName.c_str(), ios::app);
    	stressTensorFile << "xx yy zz xy yz zx \n";
    	stressTensorFile.close();

		normalStressContFileName = workDirectory + "/normalStressCont.dat";
    	normalStressContFile.open(normalStressContFileName.c_str(), ios::app);
    	normalStressContFile << "hydro-xx hydro-yy hydro-zz lub-xx lub-yy lub-zz contact-xx contact-yy contact-zz rep-xx rep-yy rep-zz \n";
    	normalStressContFile.close();

		shearStressContFileName = workDirectory + "/shearStressCont.dat";
    	shearStressContFile.open(shearStressContFileName.c_str(), ios::app);
    	shearStressContFile << "hydro-xy hydro-yz hydro-zx lub-xy lub-yz lub-zx contact-xy contact-yz contact-zx rep-xy rep-yz rep-zx \n";
    	shearStressContFile.close();

	}

    lastScreenExp = 0;
    lastFluidExp = 0;
	lastStressExp = 0;
	lastDispExp = 0;
    lastFluidLagrangianExp = 0;
	lastPartExp = 0;
    lastPart_yDispExp = 0;
	lastPart_rDispExp = 0;
	lastPart_n1Exp = 0;
	lastPart_n2Exp = 0;
	lastPart_yyStressExp = 0;
	lastPart_pressureExp = 0;
	lastBoundExp = 0;
	lastChainExp = 0;
	lastVelExp = 0;
	lastObjectExp = 0;
  lastBodyExp = 0;
  lastSpringExp = 0;
  	lastObjectExp2 = 0;
    lastFluidRecycleExp = 0;
    lastPartRecycleExp = 0;
	lastSimComplexExp = 0;
}

void IO::outputStep(const LB& lb,DEM& dem) {

    // PLOTTING PHASE  ////////////////////////////////////////////////////////////////
    const unsigned int screenExpCounter = (screenExpTime > 0 ? static_cast<unsigned int> (realTime / screenExpTime) + 1 : 0);

	if (screenExpCounter > lastScreenExp) {

        lastScreenExp = screenExpCounter;

        exportFile.open(exportFileName.c_str(), ios::app);

        // current iteration and real time
        cout << currentTimeStep << "; time=" << realTime << "\t";
        exportFile << currentTimeStep << "; time=" << realTime << "\t";

		if (dem.objects.size()){

		exportObjectPosition(dem.objects, lb);

		exportObjectForce(dem.objects, lb);


		}

    if (dem.bodies.size()){

      exportBodyPosition(dem.bodies);
    }

		if (dem.walls.size()){

		exportWallForce(dem.walls);

		}

		if (demSolve){

		exportStressTensor(dem);

		exportNormalStressContributions(dem);

		exportShearStressContributions(dem);

		}

        // closing file
        cout << "\n";
        exportFile << "\n";
        exportFile.close();
        cout.flush();
    }

    // FILE CREATION PHASE  ////////////////////////////////////////////////////////////////
    createFiles(lb, dem);

}

void IO::outputFinal() {
    // drag file closure
    exportFile.close();
}

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

// paraview files

void IO::createFiles(const LB& lb, const DEM& dem) {
    // write vtk at regular interval defined with the input file

    if (lbmSolve) {

	    const unsigned int fluidExpCounter = (fluidExpTime > 0 ? static_cast<unsigned int> (realTime / fluidExpTime) + 1 : 0);
        if (fluidExpCounter > lastFluidExp) {
            lastFluidExp = fluidExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidFileFormat.c_str(), currentTimeStep);
            exportEulerianParaviewFluid(lb, filePathBuffer);
        }

        const unsigned int fluidLagrangianExpCounter = (fluidLagrangianExpTime > 0 ? static_cast<unsigned int> (realTime / fluidLagrangianExpTime) + 1 : 0);
        if (fluidLagrangianExpCounter > lastFluidLagrangianExp) {
            lastFluidLagrangianExp = fluidLagrangianExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidLagrangianFileFormat.c_str(), currentTimeStep);
            //if (currentTimeStep>120) {
            exportLagrangianParaviewFluid(lb, filePathBuffer);
            //}
        }
        const unsigned int fluidRecycleExpCounter = (fluidRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / fluidRecycleExpTime) + 1 : 0);
        if (fluidRecycleExpCounter > lastFluidRecycleExp) {
            lastFluidRecycleExp = fluidRecycleExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, fluidRecycleFileFormat.c_str(), currentTimeStep);
            // requires the pbcShift, contained int he neighborList function
            exportRecycleFluid(lb, filePathBuffer);
        }
    }

    if (demSolve) {

        const unsigned int velExpCounter = (velExpTime > 0 ? static_cast<unsigned int> (realTime / velExpTime) + 1 : 0);
        if (velExpCounter > lastVelExp) {
            lastVelExp = velExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, velFileFormat.c_str(), currentTimeStep);
            exportParaviewParticleVelocities(dem.particles, filePathBuffer);

	    }


        const unsigned int partExpCounter = (partExpTime > 0 ? static_cast<unsigned int> (realTime / partExpTime) + 1 : 0);
        if (partExpCounter > lastPartExp) {
            lastPartExp = partExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, partFileFormat.c_str(), currentTimeStep);
            exportParaviewParticles(dem.particles, filePathBuffer);

	    }

        const unsigned int boundExpCounter = (boundExpTime > 0 ? static_cast<unsigned int> (realTime / boundExpTime) + 1 : 0);
        if (boundExpCounter > lastBoundExp) {
            lastBoundExp = boundExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, boundFileFormat.c_str(), currentTimeStep);
            exportParaviewBoundaries(dem, filePathBuffer);
        }
	if (dem.objects.size()){
        const unsigned int objectExpCounter = (objectExpTime > 0 ? static_cast<unsigned int> (realTime / objectExpTime) + 1 : 0);
        if (objectExpCounter > lastObjectExp) {
            lastObjectExp = objectExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, objectFileFormat.c_str(), currentTimeStep);
            exportParaviewObjects(dem.objects, filePathBuffer);
        }
      }

        const unsigned int bodyExpCounter = (bodyExpTime > 0 ? static_cast<unsigned int> (realTime / bodyExpTime) + 1 : 0);
        if (bodyExpCounter > lastBodyExp) {
            lastBodyExp = bodyExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, bodyFileFormat.c_str(), currentTimeStep);
            exportParaviewBodies(dem.bodies, filePathBuffer);
        }

        const unsigned int springExpCounter = (springExpTime > 0 ? static_cast<unsigned int> (realTime / springExpTime) + 1 : 0);
        if (springExpCounter > lastSpringExp) {
            lastSpringExp = springExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, springFileFormat.c_str(), currentTimeStep);
            exportParaviewSprings(dem.bodies, dem.objects, filePathBuffer);
        }


		const unsigned int partRecycleExpCounter = (partRecycleExpTime > 0 ? static_cast<unsigned int> (realTime / partRecycleExpTime) + 1 : 0);
        if (partRecycleExpCounter > lastPartRecycleExp) {
            lastPartRecycleExp = partRecycleExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, partRecycleFileFormat.c_str(), partRecycleExpCounter);
            exportRecycleParticles(dem.particles, filePathBuffer);
        }

        const unsigned int chainExpCounter = (chainExpTime > 0 ? static_cast<unsigned int> (realTime / chainExpTime) + 1 : 0);
        if (chainExpCounter > lastChainExp) {
            lastChainExp = chainExpCounter;
            char filePathBuffer [1024];
            sprintf(filePathBuffer, forceChainsFileFormat.c_str(), currentTimeStep);
            exportForceChains(dem, filePathBuffer);
        }


        }



}


void IO::exportParaviewParticles(const particleList& particles, const string& particleFile) {

    const int one = 1;
    const int Pnumber = particles.size();
    const tVect nx(1.0, 0.0, 0.0), ny(0.0, 1.0, 0.0), nz(0.0, 0.0, 1.0);
    tVect n1, n2, n3;

     std::cout.precision(10);
     std::cout.fixed;

    // file opening
    ofstream paraviewParticleFile;
    paraviewParticleFile.open(particleFile.c_str());
    // writing on header file
    paraviewParticleFile << "# vtk DataFile Version 3.0\n";
    paraviewParticleFile << "Data.vtk\n";
    paraviewParticleFile << "ASCII\n";
    paraviewParticleFile << "DATASET POLYDATA\n";
    paraviewParticleFile << "  \n";

	paraviewParticleFile << "POINTS "<< Pnumber << " float\n";
    for (int i = 0; i < Pnumber; ++i) {
        particles[i].x0.printFixedLine(paraviewParticleFile);
    }

	paraviewParticleFile << "\n";


	paraviewParticleFile << "POINT_DATA " << Pnumber << "\n";

	paraviewParticleFile << "SCALARS diameter float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << 2.0*particles[i].r << "\n";
    }

    paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS z-displacement float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].disp.y << "\n";
    }

	paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS x-displacement float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].disp.x << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS y-displacement float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].disp.z << "\n";
    }

  paraviewParticleFile << "\n";


  paraviewParticleFile << "SCALARS z-vel float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].x1.y << "\n";
    }

  paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS r-displacement float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << particles[i].disp_r << "\n";
    }

	paraviewParticleFile << "\n";


	paraviewParticleFile << "SCALARS sig-xx float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << particles[i].STotal.m00 << "\n";
    }

	paraviewParticleFile << "\n";


	paraviewParticleFile << "SCALARS sig-yy float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << particles[i].STotal.m22 << "\n";
    }

	paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS pressure float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << particles[i].pressure << "\n";
    }

	paraviewParticleFile << "\n";


	paraviewParticleFile << "SCALARS zz-stress float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].STotal.m11 << "\n";
    }

	paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS xz-stress float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].STotal.m01 + particles[i].STotal.m10)*0.5  << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS xx-stress-el float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

    for (int i = 0; i < Pnumber; ++i) {
          paraviewParticleFile << -particles[i].SElas.m00 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS yy-stress-el float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
          paraviewParticleFile << -particles[i].SElas.m22 << "\n";
    }

	paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS zz-stress-el float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].SElas.m11 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS xz-stress-el float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].SElas.m01 + particles[i].SElas.m10)*0.5  << "\n";
    }

    paraviewParticleFile << "\n";

    paraviewParticleFile << "SCALARS yz-stress-el float\n";

    paraviewParticleFile << "LOOKUP_TABLE default \n";

    for (int i = 0; i < Pnumber; ++i) {
          paraviewParticleFile << - (particles[i].SElas.m12 + particles[i].SElas.m21)*0.5  << "\n";
      }

      paraviewParticleFile << "\n";

      paraviewParticleFile << "SCALARS xy-stress-el float\n";

      paraviewParticleFile << "LOOKUP_TABLE default \n";

      for (int i = 0; i < Pnumber; ++i) {
            paraviewParticleFile << - (particles[i].SElas.m02 + particles[i].SElas.m20)*0.5  << "\n";
        }


  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS xx-stress-vis float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].SVisc.m00 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS yy-stress-vis float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].SVisc.m22 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS zz-stress-vis float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].SVisc.m11 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS xz-stress-vis float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].SVisc.m01 + particles[i].SVisc.m10)*0.5  << "\n";
    }
    paraviewParticleFile << "\n";

    paraviewParticleFile << "SCALARS xy-stress-vis float\n";

    paraviewParticleFile << "LOOKUP_TABLE default \n";

    for (int i = 0; i < Pnumber; ++i) {
          paraviewParticleFile << - (particles[i].SVisc.m02 + particles[i].SVisc.m20)*0.5  << "\n";
      }

      paraviewParticleFile << "\n";

      paraviewParticleFile << "SCALARS yz-stress-vis float\n";

      paraviewParticleFile << "LOOKUP_TABLE default \n";

      for (int i = 0; i < Pnumber; ++i) {
            paraviewParticleFile << - (particles[i].SVisc.m12 + particles[i].SVisc.m21)*0.5  << "\n";
        }


  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS zz-stress-vis-h float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << -particles[i].SVisc_h.m11 << "\n";
    }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS xz-stress-vis-h float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].SVisc_h.m01 + particles[i].SVisc_h.m10)*0.5  << "\n";
    }

  paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS xy-stress float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].STotal.m02 + particles[i].STotal.m20)*0.5  << "\n";
    }

	paraviewParticleFile << "\n";

	paraviewParticleFile << "SCALARS yz-stress float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < Pnumber; ++i) {
        paraviewParticleFile << - (particles[i].STotal.m12 + particles[i].STotal.m21)*0.5  << "\n";
    }

	paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS delta-z float\n";

	paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
    paraviewParticleFile << particles[i].overlap.y << "\n";
  }

	paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS delta-x float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
    paraviewParticleFile << particles[i].overlap.x << "\n";
  }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS delta-y float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
    paraviewParticleFile << particles[i].overlap.z << "\n";
  }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS coord float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
    paraviewParticleFile << particles[i].coord << "\n";
  }

  paraviewParticleFile << "\n";

  paraviewParticleFile << "SCALARS lub_coord float\n";

  paraviewParticleFile << "LOOKUP_TABLE default \n";

  for (int i = 0; i < Pnumber; ++i) {
    paraviewParticleFile << particles[i].lub_coord << "\n";
  }

  paraviewParticleFile << "\n";

	paraviewParticleFile.close();

}

void IO::exportParaviewParticleVelocities(const particleList& particles, const string& velFile) {

    const int one = 1;
    const int Pnumber = particles.size();
    const tVect nx(1.0, 0.0, 0.0), ny(0.0, 1.0, 0.0), nz(0.0, 0.0, 1.0);
    tVect n1, n2, n3;

     std::cout.precision(10);
     std::cout.fixed;

    // file opening
    ofstream paraviewVelFile;
    paraviewVelFile.open(velFile.c_str());
    // writing on header file
    paraviewVelFile << "# vtk DataFile Version 3.0\n";
    paraviewVelFile << "Data.vtk\n";
    paraviewVelFile << "ASCII\n";
    paraviewVelFile << "DATASET POLYDATA\n";
 //   paraviewVelFile << "  \n";

	paraviewVelFile << "POINTS "<< Pnumber << " float\n";
    for (int i = 0; i < Pnumber; ++i) {
        particles[i].x0.printFixedLine(paraviewVelFile);
    }

	paraviewVelFile << "\n";

//	paraviewVelFile << "CELL_TYPES " << Pnumber << "\n";

//	for (int i = 0; i < Pnumber ; ++i){
//		paraviewVelFile << "1\n";
//	}

//	paraviewVelFile << "\n";

	paraviewVelFile << " POINT_DATA " << Pnumber << "\n";

	paraviewVelFile << "VECTORS velocity float\n";

//	paraviewVelFile << "LOOKUP_TABLE default \n";

    for (int i = 0; i < Pnumber; ++i) {
        particles[i].x1.printFixedLine(paraviewVelFile);
    }

	paraviewVelFile << "\n";

	paraviewVelFile.close();

}

void IO::exportParaviewObjects(const objectList& objects, const string& objectFile) {

  const int one = 1;
  const int Onumber = objects.size();
  const tVect nx(1.0, 0.0, 0.0), ny(0.0, 1.0, 0.0), nz(0.0, 0.0, 1.0);
  tVect n1, n2, n3;

  std::cout.precision(10);
  std::cout.fixed;

  // file opening
  ofstream paraviewObjectFile;
  paraviewObjectFile.open(objectFile.c_str());
  // writing on header file
  paraviewObjectFile << "# vtk DataFile Version 3.0\n";
  paraviewObjectFile << "Data.vtk\n";
  paraviewObjectFile << "ASCII\n";
  paraviewObjectFile << "DATASET UNSTRUCTURED_GRID\n";
  paraviewObjectFile << "  \n";

paraviewObjectFile << "POINTS "<< Onumber << " float\n";
  for (int i = 0; i < Onumber; ++i) {
      objects[i].x0.printFixedLine(paraviewObjectFile);
  }

paraviewObjectFile << "\n";

paraviewObjectFile << "CELL_TYPES " << Onumber << "\n";

for (int i = 0; i < Onumber ; ++i){
  paraviewObjectFile << "1\n";
}

paraviewObjectFile << "\n";

paraviewObjectFile << "POINT_DATA " << Onumber << "\n";

paraviewObjectFile << "SCALARS diameter float\n";

paraviewObjectFile << "LOOKUP_TABLE default \n";

for (int i = 0; i < Onumber; ++i) {
      paraviewObjectFile << 2.0*objects[i].r << "\n";
  }

paraviewObjectFile << "\n";

paraviewObjectFile.close();

}

void IO::exportParaviewBodies(const bodyList& bodies, const string& bodyFile){

  ofstream paraviewBodyFile;
  paraviewBodyFile.open(bodyFile.c_str());
  // writing on header file
  paraviewBodyFile << "# vtk DataFile Version 3.0\n";
  paraviewBodyFile << "Data.vtk\n";
  paraviewBodyFile << "ASCII\n";
  paraviewBodyFile << "DATASET POLYDATA\n";
  paraviewBodyFile << "\n";

  paraviewBodyFile << "POINTS "<< bodies.size() << " float\n";
    for (int i = 0; i < bodies.size(); ++i) {
        bodies[i].x0.printFixedLine(paraviewBodyFile);
    }

	paraviewBodyFile << "\n";

	paraviewBodyFile << "POINT_DATA " << bodies.size() << "\n";

	paraviewBodyFile << "SCALARS diameter float\n";

	paraviewBodyFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < bodies.size(); ++i) {
        paraviewBodyFile << 2.0*bodies[i].r << "\n";
    }

    paraviewBodyFile << "\n";
    paraviewBodyFile.close();

}

void IO::exportParaviewSprings(const bodyList& bodies, const objectList& objects, const string& springFile){

  ofstream paraviewSpringFile;
  paraviewSpringFile.open(springFile.c_str());
  // writing on header file
  paraviewSpringFile << "# vtk DataFile Version 3.0\n";
  paraviewSpringFile << "Data.vtk\n";
  paraviewSpringFile << "ASCII\n";
  paraviewSpringFile << "DATASET POLYDATA\n";
  paraviewSpringFile << "\n";

  paraviewSpringFile << "POINTS "<< 2 << " float\n";
  for (int i = 0; i < bodies.size(); ++i) {
      bodies[i].x0.printFixedLine(paraviewSpringFile);
  }
  for (int i = 0; i < objects.size(); ++i) {
      objects[i].x0.printFixedLine(paraviewSpringFile);
  }
  paraviewSpringFile << "\n";

  paraviewSpringFile << "LINES " << 1 << " " << 3 << "\n";
  paraviewSpringFile <<"2 " << 0 << " " << 1 << "\n";

  paraviewSpringFile << "\n";

  paraviewSpringFile << "POINT_DATA " << 2 << "\n";
  paraviewSpringFile << "SCALARS radius float\n";
  paraviewSpringFile << "LOOKUP_TABLE default \n";
  paraviewSpringFile << "1.0\n";
  paraviewSpringFile << "1.0\n";
  paraviewSpringFile << "\n";
  paraviewSpringFile << "SCALARS compression float\n";
  paraviewSpringFile << "LOOKUP_TABLE default \n";
  for (int i = 0; i < bodies.size(); ++i) {
        paraviewSpringFile << bodies[i].L - bodies[i].L0 << "\n";
        paraviewSpringFile << bodies[i].L - bodies[i].L0 << "\n";
    }
  paraviewSpringFile << "\n";
  paraviewSpringFile.close();
}



void IO::exportParaviewBoundaries(const DEM& dem,  const string& boundaryFile) {

     std::cout.precision(10);
     std::cout.fixed;


    // file opening
    ofstream paraviewBoundaryFile;
    paraviewBoundaryFile.open(boundaryFile.c_str());

	// writing on header file
    paraviewBoundaryFile << "# vtk DataFile Version 3.0\n";
    paraviewBoundaryFile << "DataFC.vtk\n";
    paraviewBoundaryFile << "ASCII\n";
    paraviewBoundaryFile << "DATASET POLYDATA\n";
    paraviewBoundaryFile << "\n";

	paraviewBoundaryFile << "POINTS "<< 24 << " float\n";

// position of corners

    paraviewBoundaryFile << 0.0  << " " <<  0.0<< " " << 0.0 << "\n";
    paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << 0.0 << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << 0.0 << "\n";
	paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << dem.demSize[2] << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  0.0 << " " << dem.demSize[2] << "\n";

	paraviewBoundaryFile << 0.0 << " " <<  0.0 << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  0.0 << " " << 0.0 << "\n";

//

	paraviewBoundaryFile << 0.0 << " " <<  0.0 << " " << 0.0 << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << 0.0 << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << 0.0 << "\n";
	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << 0.0 << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  0.0 << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";

	paraviewBoundaryFile << 0.0 << " " <<  0.0 << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";

//
	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << 0.0 << "\n";
	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << 0.0 << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << 0.0 << "\n";
	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";

	paraviewBoundaryFile << dem.demSize[0] << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";

	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << dem.demSize[2] << "\n";
	paraviewBoundaryFile << 0.0 << " " <<  dem.demSize[1] << " " << 0.0 << "\n";


	paraviewBoundaryFile << "\n";

	paraviewBoundaryFile << "LINES " << 12 << " " << 36 << "\n";

	for (int i = 0 ; i < 12 ; i++){
		paraviewBoundaryFile <<"2 " << 2*i << " " << 2*i + 1 << "\n";
	}

	paraviewBoundaryFile << "\n";

	paraviewBoundaryFile << "POINT_DATA " << 24 << "\n";

	paraviewBoundaryFile << "SCALARS radius float\n";

	paraviewBoundaryFile << "LOOKUP_TABLE default \n";

	for (int i = 0; i < 24 ; ++i){
		paraviewBoundaryFile << "1.0\n";
	}

	paraviewBoundaryFile.close();

}

void IO::exportForceChains(const DEM& dem,  const string& forceChainsFile) {

    std::cout.precision(10);
    std::cout.fixed;

	int n_con = dem.numContact;

    // file opening
    ofstream paraviewForceChainsFile;
    paraviewForceChainsFile.open(forceChainsFile.c_str());

	// writing on header file
    paraviewForceChainsFile << "# vtk DataFile Version 3.0\n";
    paraviewForceChainsFile << "DataFC.vtk\n";
    paraviewForceChainsFile << "ASCII\n";
    paraviewForceChainsFile << "DATASET POLYDATA\n";
    paraviewForceChainsFile << "\n";

	paraviewForceChainsFile << "POINTS "<< 2*(n_con+1) << " float\n";
	paraviewForceChainsFile << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	paraviewForceChainsFile << 0.0 << " " << 0.0 << " " << 0.0 << endl;

	for (int n = 0 ;  n < n_con; n++){
		paraviewForceChainsFile << dem.cPointIx[n] << " " << dem.cPointIy[n] << " " << dem.cPointIz[n] << endl;
		paraviewForceChainsFile << dem.cPointJx[n] << " " << dem.cPointJy[n] << " " << dem.cPointJz[n] << endl;


	}

	paraviewForceChainsFile << "\n";

	paraviewForceChainsFile << "LINES " << n_con+1 << " " << 3*(n_con+1)<< "\n";

// position of corners
	for (int i = 0 ; i < n_con+1 ; i++){
	 paraviewForceChainsFile <<"2 " << 2*i << " " << 2*i + 1 << "\n";
	}

	paraviewForceChainsFile << "\n";


	paraviewForceChainsFile << "POINT_DATA " << 2*(n_con+1) << "\n";

	paraviewForceChainsFile << "SCALARS radius float\n";

	paraviewForceChainsFile << "LOOKUP_TABLE default \n";

	paraviewForceChainsFile << 0.0 << endl;

	paraviewForceChainsFile << 0.2*dem.sphereMat.linearStiff << endl;

	for (int i = 0; i < n_con ; ++i){
		paraviewForceChainsFile << dem.contactForceN[i] << endl;
		paraviewForceChainsFile << dem.contactForceN[i] << endl;
	}
	paraviewForceChainsFile.close();
}


void IO::exportRecycleFluid(const LB& lb, const string& fluidRecycleFile) {
    // exports a format readable by this code itself, allowing the re-use of computed particle data.
    // needs further work...

    ofstream recycleFluidFile;
    recycleFluidFile.open(fluidRecycleFile.c_str());

    recycleFluidFile << lb.lbSize[0] << " " << lb.lbSize[1] << " " << lb.lbSize[2] << "\n";

    for (int it = 0; it < lb. totNodes; ++it) {
        recycleFluidFile << lb.types[it].getType() << " ";
    }
    recycleFluidFile << endl;
    cout << endl;
    for (int it = 0; it < lb.activeNodes.size(); ++it) {
        const int index = lb.activeNodes[it];
        // to gain precision, scale density and mass by 1 (initDensity)
        recycleFluidFile << lb.nodes[index]->n - lb.fluidMaterial.initDensity << " " << lb.nodes[index]->u.dot(Xp) << " " << lb.nodes[index]->u.dot(Yp) << " " << lb.nodes[index]->u.dot(Zp) << " " << lb.nodes[index]->mass - lb.fluidMaterial.initDensity << " " << lb.nodes[index]->visc << " ";
        for (int j = 0; j < lbmDirec; ++j) {
            recycleFluidFile << lb.nodes[index]->f[j] << " ";
        }
        recycleFluidFile << endl;
    }
    recycleFluidFile.close();
}

void IO::exportRecycleParticles(const particleList& particles,  const string& partRecycleFile) {
	std::cout.precision(10);

    ofstream recycleParticleFile;
    recycleParticleFile.open(partRecycleFile.c_str());

   // recycleParticleFile << particles.size() << "\n";

    for (int i = 0; i < particles.size(); i++) {
        // import variables
        recycleParticleFile << particles[i].r << " ";

        particles[i].x0.print(recycleParticleFile);

		particles[i].x1.print(recycleParticleFile);

		recycleParticleFile << -particles[i].STotal.m00 << "\t";
		recycleParticleFile << -particles[i].STotal.m11 << "\t";
		recycleParticleFile << -particles[i].STotal.m22 << "\t";
		recycleParticleFile << -(particles[i].STotal.m01 + particles[i].STotal.m10)*0.5 << "\t";
		recycleParticleFile << -(particles[i].STotal.m12 + particles[i].STotal.m21)*0.5 << "\t";
		recycleParticleFile << -(particles[i].STotal.m02 + particles[i].STotal.m20)*0.5 << "\t";
		recycleParticleFile << -particles[i].disp.y  << "\n";


      //  particles[i].x1.print(recycleParticleFile);

	/*	recycleParticleFile << particles[i].FHydro.norm() << " ";
		recycleParticleFile << particles[i].FLub.norm() << " ";
		recycleParticleFile << particles[i].FContact.norm() << " ";
		recycleParticleFile << particles[i].FGrav.norm() << " ";
		recycleParticleFile << particles[i].FRep.norm() << " ";
		recycleParticleFile << particles[i].FTotal.norm() << " ";
		recycleParticleFile << " \n"; */

    }
    recycleParticleFile.close();
}

void IO::exportLagrangianParaviewFluid(const LB& lb, const string& fluidFile) {

    const int one = 1;

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;
    paraviewFluidFile.open(fluidFile.c_str());
    // writing on header file
    paraviewFluidFile << "<?xml version=\"1.0\"?>\n";
    paraviewFluidFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    paraviewFluidFile << " <UnstructuredGrid GhostLevel=\"0\">\n";
    paraviewFluidFile << "  <Piece NumberOfPoints=\"" << lb.activeNodes.size() << "\" NumberOfCells=\"" << lb.activeNodes.size() << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
        tVect uPhysic = lb.nodes[index]->u * lb.unit.Speed;
        uPhysic.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
            paraviewFluidFile << 0.3333333 * (lb.nodes[index]->n - lb.fluidMaterial.initDensity) * lb.unit.Pressure << "\n";
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"n\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
            paraviewFluidFile << lb.nodes[index]->n << "\n";
    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"visc\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
            paraviewFluidFile << lb.nodes[index]->visc << "\n";
    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"viscosity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (i = 0; i < Pnumber; ++i) {
    //        it = lb.activeNodes[i];
    //        paraviewFluidFile << lb.nodes[it]->visc * lb.unit.KinVisc << "\n";
    //    }
    //    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"isInit\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (i=0; i<Pnumber; ++i){
    //
    //        unsigned int a=0;
    //        if (lb.nodes[i]!=0) {
    //            a=1;
    //        }
    //        paraviewFluidFile<<a<<"\n";
    //    }
    //    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"shearRate\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    //    for (i = 0; i < Pnumber; ++i) {
    //        it = lb.activeNodes[i];
    //        paraviewFluidFile << lb.nodes[it]->shearRate * lb.unit.AngVel << "\n";
    //    }
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
        paraviewFluidFile << lb.nodes[index]->mass * lb.unit.Density << "\n";
    }
    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
        paraviewFluidFile << lb.types[index].getType() << "\n";
    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "   <Points>\n";
    paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        const unsigned int index = lb.activeNodes[i];
        const tVect positionHere = lb.getPosition(index)* lb.unit.Length;
        positionHere.printLine(paraviewFluidFile);
    }
    paraviewFluidFile << "   </Points>\n";
    paraviewFluidFile << "   <Cells>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 1; i < lb.activeNodes.size() + 1; ++i) {
        paraviewFluidFile << i - 1 << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 1; i < lb.activeNodes.size() + 1; ++i) {
        paraviewFluidFile << i << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "    <DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < lb.activeNodes.size(); ++i) {
        paraviewFluidFile << one << "\n";
    }
    paraviewFluidFile << "    </DataArray>\n";
    paraviewFluidFile << "   </Cells>\n";
    paraviewFluidFile << "  </Piece>\n";
    //    paraviewHeaderFile<<"  <Piece Source=\""<<pieceFile<<"\"/>\n";
    paraviewFluidFile << " </UnstructuredGrid>\n";
    paraviewFluidFile << "</VTKFile>";

    // data file closing
    paraviewFluidFile.close();
}

void IO::exportEulerianParaviewFluid(const LB& lb, const string& fluidFile) {

    //    const char* charFluidFile;
    //
    //    string title, fluidFile;
    //    title=fluidDirectory+"/fluid_";
    //    stringstream ss;
    //    ss<<currentTimeStep;
    //    fluidFile=title+ss.str()+".vti";
    //    charFluidFile=fluidFile.c_str();

    const double zero = 0.0;

    // start printing all the crap required for Paraview
    // header file opening
    ofstream paraviewFluidFile;
    paraviewFluidFile.open(fluidFile.c_str());
    // writing on header file
    paraviewFluidFile << "<VTKFile type=\"ImageData\" version=\"0.1\">\n";
    paraviewFluidFile << " <ImageData WholeExtent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\" "
            << "Origin=\"0.0 0.0 0.0\" Spacing=\"" << lb.unit.Length << " " << lb.unit.Length << " " << lb.unit.Length << "\">\n";
    paraviewFluidFile << "  <Piece Extent=\"0 " << lb.lbSize[0] - 1 << " 0 " << lb.lbSize[1] - 1 << " 0 " << lb.lbSize[2] - 1 << "\">\n";
    paraviewFluidFile << "   <PointData>\n";
    paraviewFluidFile << "    <DataArray type=\"Int8\" Name=\"type\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"8\">\n";
    for (int i = 0; i < lb.totNodes; ++i) {
	       if (!lb.types[i].isInsideParticle()) {
			           paraviewFluidFile << lb.types[i].getType() << " ";
				      } else {
		          paraviewFluidFile << 1 << " ";
				       }
			   }
    paraviewFluidFile << "    </DataArray>\n";
    //    for (int j=1; j<lbmDirec; ++j) {
    //        paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"j"<<j<<"\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    //        for (i=0; i<Pnumber; ++i){
    //            if (lb.curves[i]==0) {
    //                Zero.print(paraviewFluidFile);
    //            }
    //            else {
    //                tVect deltaVec=lb.unit.Length*lb.curves[i]->delta[j]*v[j];
    //                deltaVec.print(paraviewFluidFile);
    //            }
    //        }
    //        paraviewFluidFile<<"    </DataArray>\n";
    //    }
 /*   paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"v\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < lb.totNodes; ++i) {
        if (lb.nodes[i] == 0) {
            Zero.print(paraviewFluidFile);
        } else {
            const tVect uPhysic = lb.nodes[i]->u * lb.unit.Speed;
            uPhysic.print(paraviewFluidFile);
        }
    }
    paraviewFluidFile << "    </DataArray>\n"; */
    //    paraviewFluidFile<<"    <DataArray type=\"Float64\" Name=\"f\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    //    for (i=0; i<Pnumber; ++i){
    //        if (lb.nodes[i]==0) {
    //            Zero.print(paraviewFluidFile);
    //        }
    //        else {
    //            tVect fPhysic=lb.nodes[i]->hydroForce*lb.unit.Force;
    //            fPhysic.print(paraviewFluidFile);
    //        }
    //    }
    //    paraviewFluidFile<<"    </DataArray>\n";
//    if (lb.freeSurface) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"amass\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totNodes; ++i) {
            if (lb.nodes[i] == 0) {
                paraviewFluidFile << zero << " ";
            } else {
                paraviewFluidFile << lb.nodes[i]->mass * lb.unit.Density << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
		// }
    if (lb.lbTopography) {
        paraviewFluidFile << "    <DataArray type=\"Float64\" Name=\"topSurface\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n";
        for (int i = 0; i < lb.totNodes; ++i) {
            if (lb.types[i].isTopography()) {
                paraviewFluidFile << zero << " ";
            } else {
                paraviewFluidFile << 1.0 << " ";
            }
        }
        paraviewFluidFile << "    </DataArray>\n";
    }
    paraviewFluidFile << "   </PointData>\n";
    paraviewFluidFile << "   <CellData>\n";
    paraviewFluidFile << "   </CellData>\n";
    paraviewFluidFile << "  </Piece>\n";
    paraviewFluidFile << " </ImageData>\n";
    paraviewFluidFile << "</VTKFile>\n";
    // data file closing
    paraviewFluidFile.close();
}

// quantities

void IO::exportObjectForce(const objectList& objects,const LB& lb){

  objectForceFile.open(objectForceFileName.c_str(), ios::app);

	for (int o = 0; o < objects.size(); o++ )
	{
		double depth = (lb.fluidHeight + objects[o].r) - objects[o].x0.y;

		tVect FExerted = tVect(0.0,0.0,0.0);

		FExerted = objects[o].FTotal - objects[o].FGrav;

		tVect SExerted = tVect(0.0,0.0,0.0);

		if (depth > 0){

		double area = 4.0 * M_PI * objects[o].r * objects[o].r * ( depth / (2.0 * objects[o].r) );

		SExerted = FExerted / area ;

		}

		FExerted.print(objectForceFile);

		objectForceFile  << "\t";

		SExerted.print(objectForceFile);

		objectForceFile  << "\t";

		objects[o].FTotal.print(objectForceFile);

		objectForceFile  << "\t";

		objects[o].FContact.print(objectForceFile);

		objectForceFile  << "\t";

		objects[o].FHydro.print(objectForceFile);

		objectForceFile  << "\t";

		objects[o].FLub.print(objectForceFile);

		objectForceFile  << "\t";

		objects[o].FRep.print(objectForceFile);

		objectForceFile << endl;

		}
	objectForceFile.close();


}

void IO::exportWallForce(const wallList& walls){

	wallForceFile.open(wallForceFileName.c_str(), ios::app);

	for (int o = 0; o < walls.size(); o++ )
	{
		walls[o].FTotal.print(wallForceFile);

		wallForceFile  << "\t";

		walls[o].FContact.print(wallForceFile);

		wallForceFile  << "\t";

		}

	wallForceFile << endl;
	wallForceFile.close();

}

void IO::exportObjectPosition(const objectList& objects, const LB& lb){

  objectPositionFile.open(objectPositionFileName.c_str(), ios::app);

	for (int o = 0; o < objects.size(); o++ )
	{

	double depth = (lb.fluidHeight + objects[o].r) - objects[o].x0.y;

	objects[o].x0.print(objectPositionFile);

	objectPositionFile  << "\t";

	objects[o].x1.print(objectPositionFile);

	objectPositionFile  << "\t" << depth;

	objectPositionFile << endl;

	}
	objectPositionFile.close();

}

void IO::exportBodyPosition(const bodyList& bodies){

  bodyPositionFile.open(bodyPositionFileName.c_str(), ios::app);
  for (int o = 0; o < bodies.size() ; o++ ){
    bodies[o].x0.print(bodyPositionFile);
    bodyPositionFile  << "\t";
    bodies[o].x1.print(bodyPositionFile);
    bodyPositionFile  << "\t" << bodies[o].L;
    bodyPositionFile << endl;
  }

  bodyPositionFile.close();

}

void IO::exportStressTensor(const DEM& dem){
	stressTensorFile.open(stressTensorFileName.c_str(), ios::app);

	dem.totalStress.print(stressTensorFile);

	stressTensorFile << endl;

	stressTensorFile.close();
}

void IO::exportNormalStressContributions(const DEM& dem){
	normalStressContFile.open(normalStressContFileName.c_str(), ios::app);

	dem.hydroStress.printNormal(normalStressContFile);

	normalStressContFile  << "\t";

	dem.lubStress.printNormal(normalStressContFile);

	normalStressContFile  << "\t";

	dem.contactStress.printNormal(normalStressContFile);

	normalStressContFile  << "\t";

	dem.electroStress.printNormal(normalStressContFile);

	normalStressContFile  << endl;

	normalStressContFile.close();

}

void IO::exportShearStressContributions(const DEM& dem){
	shearStressContFile.open(shearStressContFileName.c_str(), ios::app);

	dem.hydroStress.printShear(shearStressContFile);

	shearStressContFile  << "\t";

	dem.lubStress.printShear(shearStressContFile);

	shearStressContFile  << "\t";

	dem.contactStress.printShear(shearStressContFile);

	shearStressContFile  << "\t";

	dem.electroStress.printShear(shearStressContFile);

	shearStressContFile  << endl;

	shearStressContFile.close();

}
