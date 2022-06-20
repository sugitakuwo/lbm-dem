
#include "DEM.h"

#define USE_MATH_DEFINES
#define  n_int(x)       ( (x) >= (0) ? (int)((x)+0.5) : (int)((x)-0.5))
#define  s_int(x)       ( (x) >= (0) ? (int)(x) : (int)((x)-1.0))
#define  sgn(x)         ( (x) <  (0) ? (-1) : (1))
#define  sign(x)        (((x) < 0 ) ? -1 : ( (x) > 0 ))
#define  max(a, b)      ( (a) >  (b) ? (a)  : (b))
#define  min(a, b)      ( (a) <  (b) ? (a)  : (b))
#define  box(a, b)      ( (a) - s_int((a)/(b))*(b))
#define  n_image(a, b)  ( (a) - n_int((a)/(b))*(b))
#define  range(x, a, b) (((x) >= (a) && (x) < (b)) ? (1) : (0))

using namespace std;

void DEM::discreteElementGet(GetPot& configFile, GetPot& commandLine) {

    PARSE_CLASS_MEMBER(configFile, sphereMat.density, "density", 0.0);
    ASSERT(sphereMat.density > 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.objectDensity, "objectDensity", 0.0);
    ASSERT(sphereMat.objectDensity > 0.0);


    PARSE_CLASS_MEMBER(configFile, sphereMat.bodyDensity, "bodyDensity", 0.0);
    ASSERT(sphereMat.bodyDensity > 0.0);


    PARSE_CLASS_MEMBER(configFile, sphereMat.linearStiff, "linearStiff", 0.0);
    ASSERT(sphereMat.linearStiff >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.frictionCoefPart, "frictionCoefPart", 0.0);
    ASSERT(sphereMat.frictionCoefPart >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.lubNormal, "lubNormal", 0.0);
    ASSERT(sphereMat.lubNormal >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.lubTangential, "lubTangential", 0.0);
    ASSERT(sphereMat.lubTangential >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.lubRotation, "lubRotation", 0.0);
    ASSERT(sphereMat.lubRotation >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.lubCutoff, "lubCutoff", 0.0);
    ASSERT(sphereMat.lubCutoff >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.kinVisc, "kinVisc", 0.0);
    ASSERT(sphereMat.kinVisc >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.electroForce, "electroForce", 0.0);
    ASSERT(sphereMat.electroForce >= 0.0);

    PARSE_CLASS_MEMBER(configFile, sphereMat.debyeLength, "debyeLength", 0.0);
    ASSERT(sphereMat.debyeLength >= 0.0);

    PARSE_CLASS_MEMBER(configFile, objWaitTime, "objWaitTime", 0);
    ASSERT(objWaitTime >= 0);

    PARSE_CLASS_MEMBER(configFile, intervalTime, "intervalTime", 0);
    ASSERT(intervalTime >= 0);

    PARSE_CLASS_MEMBER(configFile, springDamp, "springDamp", 0.0);
    ASSERT(springDamp >= 0.0);

    PARSE_CLASS_MEMBER(configFile, epsilon, "epsilon", 0.0);
    ASSERT(epsilon >= 0.0);

    PARSE_CLASS_MEMBER(configFile, nonLinType, "nonLinType", 0);
    ASSERT(nonLinType >= 0);

    PARSE_CLASS_MEMBER(configFile, k13, "k13", 0.0);
    ASSERT(k13 >= 0.0);



    double scale = 1.0;
    PARSE_CLASS_MEMBER(configFile, scale, "scale", 1.0);

    ifstream particleFileID;
    particleFileID.open("./input/particles.dat");
    cout<<"Reading "<<"./input/particles.dat"<<endl;
    ASSERT(particleFileID.is_open());

    unsigned int totElmt;
    particleFileID>>totElmt;

    for (int n = 0; n < totElmt; ++n) {
        particle dummyElmt;

        // import variables
        particleFileID >> dummyElmt.index;
        particleFileID >> dummyElmt.size;
        particleFileID >> dummyElmt.r;
        dummyElmt.r = dummyElmt.r*scale;
        // position
        double x0, y0, z0;
        particleFileID>>x0;
        particleFileID>>y0;
        particleFileID>>z0;
        dummyElmt.x0 = tVect(x0*scale - 0.05, y0*scale - 0.05, z0 * scale - 0.05) ;
        // translational velocity
        double x1, y1, z1;
        particleFileID>>x1;
        particleFileID>>y1;
        particleFileID>>z1;
        dummyElmt.x1 = tVect(x1, y1, z1); //.reset();
        // rotational velocity
        double w01, w02, w03;
        particleFileID>>w01;
        particleFileID>>w02;
        particleFileID>>w03;
        dummyElmt.w1 = tVect(w01, w02, w03); //.reset();
        //dummyElmt.x1.reset();
        // orientation
        double p0, q0, r0, s0;
        particleFileID>>p0;
        particleFileID>>q0;
        particleFileID>>r0;
        particleFileID>>s0;
        dummyElmt.q0 = tQuat(p0, q0, r0, s0);
        //dummyElmt.q0.resetSoft();
        // translational velocity (in quaternion rates))
        double p1, q1, r1, s1;
        particleFileID>>p1;
        particleFileID>>q1;
        particleFileID>>r1;
        particleFileID>>s1;
        dummyElmt.q1 = tQuat(p1, q1, r1, s1);
        //dummyElmt.q1.resetHard();
        // add to list
        particles.push_back(dummyElmt);

    }


    ifstream objectFileID;
    objectFileID.open("./input/objects.dat");
    cout<<"Reading "<<"./input/objects.dat"<<endl;
    ASSERT(objectFileID.is_open());

    unsigned int totObject;
    objectFileID>>totObject;

    for (int n = 0; n < totObject; ++n) {
        object dummyObject;

        // import variables
        objectFileID >> dummyObject.index;
        objectFileID >> dummyObject.r;
        dummyObject.r = dummyObject.r*scale;
        // position
        double x0, y0, z0;
        objectFileID>>x0;
        objectFileID>>y0;
        objectFileID>>z0;
        dummyObject.x0 = tVect(x0*scale, y0*scale, z0 * scale) ;
        // translational velocity
        double x1, y1, z1;
        objectFileID>>x1;
        objectFileID>>y1;
        objectFileID>>z1;
        dummyObject.x1 = tVect(x1, y1, z1); //.reset();
		objects.push_back(dummyObject);

	}

  ifstream bodyFileID;
  bodyFileID.open("./input/body.dat");
  cout<<"Reading "<<"./input/body.dat"<<endl;
  ASSERT(bodyFileID.is_open());

  unsigned int totBody;
  bodyFileID>>totBody;

  for (int n = 0; n < totBody; ++n) {
      body dummyBody;

      // import variables
      bodyFileID >> dummyBody.index;
      bodyFileID >> dummyBody.r;
      double x0, y0, z0;
      bodyFileID>>x0;
      bodyFileID>>y0;
      bodyFileID>>z0;
      dummyBody.x0 = tVect(x0*scale, y0*scale, z0 * scale) ;
      // translational velocity
      double x1, y1, z1;
      bodyFileID>>x1;
      bodyFileID>>y1;
      bodyFileID>>z1;
      dummyBody.x1 = tVect(x1, y1, z1); //.reset();
      bodyFileID >> dummyBody.k;
      bodyFileID >> dummyBody.L0;
      bodies.push_back(dummyBody);
    }
}

void DEM::discreteElementInit(const unsIntList& externalBoundary, const doubleList& externalSize, const vecList& externalBoundaryLocation,  const double& externalTimeStep, const tVect externalAccel, const double& height, bool& surface) {

    demSize[0] = externalSize[0];

	if (surface){
		demSize[1] = height;

	}else{
		demSize[1] = externalSize[1];
	}

    demSize[2] = externalSize[2];

	if (surface){
		vol = demSize[0]*height*demSize[2];
	}else{
		vol = demSize[0]*demSize[1]*demSize[2];
	}

    // initializing particles
    const double partDensity = sphereMat.density;
	  const double objDensity = sphereMat.objectDensity;
    const double bodDensity = sphereMat.bodyDensity;
    // acceleration field
    demF = externalAccel;

	unsigned int globalIndex = 0;



 	tVect zeroGrav = tVect(0.0,0.0,0.0);

	for (int n = 0; n < particles.size(); n++) {
        particles[n].initialize(partDensity, demF, globalIndex);
    }

	// initialize objects

  for (int n = 0; n < objects.size(); n++) {


		objects[n].m = 4.0/3.0*objDensity*M_PI*objects[n].r*objects[n].r*objects[n].r;

		objects[n].FGrav = demF * objects[n].m;

		objects[n].I = 2.0/5.0 * objects[n].m * objects[n].r* objects[n].r;

	 	objects[n]._m = 1.0 / objects[n].m;

		objects[n]._I = 1.0 / objects[n].I;
    objects[n].rho = sphereMat.objectDensity;

	 	objects[n]._x0.x = objects[n].x0.x - deltat*objects[n].x1.x ;
	 	objects[n]._x0.y = objects[n].x0.y - deltat*objects[n].x1.y ;
	 	objects[n]._x0.z = objects[n].x0.z - deltat*objects[n].x1.z ;

	}

  for (int n = 0; n < bodies.size(); n++) {

    bodies[n].m = bodDensity*(4.0/3.0)*M_PI*(bodies[n].r*bodies[n].r*bodies[n].r);

    bodies[n]._m = 1.0 / bodies[n].m;

    bodies[n].FGrav = demF*bodies[n].m;

    bodies[n].L = bodies[n].L0;

  }


    // initializing wall for DEM
  initializeWalls(externalBoundary, externalBoundaryLocation, globalIndex);

	// initialization for verlet integration

	for (int n = 0; n < particles.size(); n++){

		 particles[n]._m = 1.0 / particles[n].m;
		 particles[n]._I = 1.0 / particles[n].I;

		 particles[n]._x0.x = particles[n].x0.x - deltat*particles[n].x1.x ;
		 particles[n]._x0.y = particles[n].x0.y - deltat*particles[n].x1.y ;
		 particles[n]._x0.z = particles[n].x0.z - deltat*particles[n].x1.z ;

		 particles[n].xt0.x = particles[n]._x0.x;
		 particles[n].xt0.y = particles[n]._x0.y;
		 particles[n].xt0.z = particles[n]._x0.z;

     particles[n].x_i = particles[n].x0;

	 }

	// initialize neighbor list parameters

	neighborList();

	newNeighborList = true ;

    double totMass(0.0);
    minPartRadius = 1.0e99;
    maxPartRadius = 0.0;
    meanPartRadius= 0.0;
    for (int n = 0; n < particles.size(); ++n) {
        const double radiusHere = particles[n].r;
        // calculate mass
        totMass += particles[n].m;
        meanPartRadius += radiusHere;
        if (radiusHere > maxPartRadius) {
            maxPartRadius = radiusHere;
        }
        if (radiusHere < minPartRadius) {
            minPartRadius = radiusHere;
        }
    }
    if (particles.size() > 0.0) {
        meanPartRadius /= double(particles.size());
    }



    cout << "DEM parameters\n";
    cout << "domain size: xdim =" << demSize[0] << "; ydim= " << demSize[1] << "; zdim= " << demSize[2] << ";"<<endl;
    cout << "Tot particles: " << particles.size() << ";\t";

    cout << "Particle radius: Mean="<<meanPartRadius<<" Min="<<minPartRadius<<" Max="<<maxPartRadius<<";"<<endl;
    cout << "Total particle mass=" << totMass << endl;

	  deltat = externalTimeStep / (double(multiStep));

    cout << "Deltat =" << deltat << ", by imposing " << multiStep << " substeps" << endl;


    cout << "Contact model: linear dashpot" << endl;
    cout << "Normal stiffness = " << sphereMat.linearStiff << endl;


    double averageMass = 4.0/3.0 * partDensity * M_PI * meanPartRadius *  meanPartRadius *  meanPartRadius ;

	  contactDuration = sqrt(averageMass/sphereMat.linearStiff);

	  sphereMat.dampCoeff = sqrt(averageMass*sphereMat.linearStiff);

    cout << "Damping ratio = " << sphereMat.dampCoeff << endl;

	  cout << "Contact duration = " << contactDuration << endl;

    cout << "Friction coefficient = " << sphereMat.frictionCoefPart << endl;

	  cout << "Magnitude of electrostatic force = " << sphereMat.electroForce << endl;

	  cout << "Debye length = " << sphereMat.debyeLength << endl;

    cout << "Waiting time for objects to move = " << objWaitTime << " "<< "timesteps" << endl;

    cout << "Interval time = " << intervalTime << " "<< "timesteps" << endl;


}

void DEM::discreteElementStep(const unsIntList& externalBoundary) {

    // set trigger for new neighbor list
    static const double neighListTrigger = 0.05 ;

    for (int demIter = 0; demIter < multiStep; ++demIter) {

        demTimeStep++;
        demTime += deltat;

        // neighbor management
        evalMaxDisp();

	    if (maxDisp > neighListTrigger)
		    {
           maxDisp = 0.0;
           cout<<"new neighbor list"<<endl;
           neighborList();
		   newNeighborList = true ;
        }

		evaluateForces(demTimeStep) ;

		integration(demTimeStep);


		if (externalBoundary[0] == 4 && externalBoundary[1] == 4){

			xPbcs();
		}

		if (externalBoundary[2] == 4 && externalBoundary[3] == 4){

			yPbcs();
		}

		if (externalBoundary[4] == 4 && externalBoundary[5] == 4){

			zPbcs();
		  }
    }

}



void DEM::initializeWalls(const unsIntList& externalBoundary, const vecList& boundaryLocation, unsigned int& globalIndex) {
    walls.clear();

    // boundary directors
    vecList boundaryDirec;
    boundaryDirec.resize(6);
    boundaryDirec[0] = tVect(1.0, 0.0, 0.0);
    boundaryDirec[1] = tVect(-1.0, 0.0, 0.0);
    boundaryDirec[2] = tVect(0.0, 1.0, 0.0);
    boundaryDirec[3] = tVect(0, 0. - 1.0, 0.0);
    boundaryDirec[4] = tVect(0.0, 0.0, 1.0);
    boundaryDirec[5] = tVect(0.0, 0.0, -1.0);

    unsigned int index = 0;
    // basic walls
    for (int i = 0; i < externalBoundary.size(); ++i) {
        if ((externalBoundary[i] == 5) || (externalBoundary[i] == 6) || (externalBoundary[i] == 7) || (externalBoundary[i] == 8)) {
            wall dummyWall;
            dummyWall.p = boundaryLocation[i]; //tVect(0.5*unit.Length,0.0,0.0);
            dummyWall.n = boundaryDirec[i];
            dummyWall.index = index;
			      dummyWall.wallIndex = globalIndex;
            dummyWall.translating = false;
            dummyWall.trans.reset();
            dummyWall.limited = false;
            if (externalBoundary[0] == 5) {
                dummyWall.moving = false;
                dummyWall.slip = true;
            } else if (externalBoundary[0] == 6) {
                dummyWall.moving = true;
                dummyWall.slip = true;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel.reset();
            } else if (externalBoundary[0] == 7) {
                dummyWall.moving = false;
                dummyWall.slip = false;
            } else if (externalBoundary[0] == 8) {
                dummyWall.moving = true;
                dummyWall.slip = false;
                dummyWall.rotCenter.reset();
                dummyWall.omega.reset();
                dummyWall.vel = tVect(0.0, 0.0, 0.0);
            }
			++globalIndex;
            ++index;
            walls.push_back(dummyWall);
        }
    }

    for (int n = 0; n < walls.size(); ++n) {
        walls[n].wallShow();
    }
}

void DEM::evalMaxDisp() {

    double maxVel = 0.0;
    for (int n = 0; n < particles.size(); n++) {
        if (particles[n].x1.norm2() > maxVel) {
            maxVel = particles[n].x1.norm2();
        }
    }
    maxVel = sqrt(maxVel);
    maxDisp += maxVel*deltat;
}

//neighbor list function, from Susp3D by Anthony Ladd
void DEM::neighborList(){

	// maximum neighbor
	const int MAX_N = 50;

	// max link cells
	const int MAX_L = 100000;

	static int cell_list[MAX_L][MAX_N+1];

	double x12, y12, z12, r12, r_max, sig, cut;
	int nlx, nly, nlz;
	int lx1, ly1, lz1, l1, n1, lx2, ly2, lz2, l2, n2;
	int n1_sph, n2_sph,n_sph;
	int n, j, l, k;

	int num_sph = particles.size();

	for (n = 0 ; n < num_sph; n++){

		memcpy(particles[n]._list,particles[n].list, (MAX_N) * sizeof(int));
		memcpy(particles[n]._xix,particles[n].xix, (MAX_N) * sizeof(double));
		memcpy(particles[n]._xiy,particles[n].xiy, (MAX_N) * sizeof(double));
		memcpy(particles[n]._xiz,particles[n].xiz, (MAX_N) * sizeof(double));
	}

	for (n1_sph = 0, r_max = 0.0; n1_sph < num_sph; n1_sph++)          /* Largest radius */
	{
		r_max = max(r_max, particles[n1_sph].r);
		particles[n1_sph].x0_lst.x = particles[n1_sph].x0.x;                 /* Store positions */
		particles[n1_sph].x0_lst.y = particles[n1_sph].x0.y;
		particles[n1_sph].x0_lst.z = particles[n1_sph].x0.z;
		particles[n1_sph].list[0] = 0;
		particles[n1_sph].xix[0] = 0.0;
  	particles[n1_sph].xiy[0] = 0.0;
		particles[n1_sph].xiz[0] = 0.0;
	}

	cut = 2.0*r_max + 2.0;
	nlx = demSize[0]/cut; nly = demSize[1]/cut; nlz = demSize[2]/cut;


	if (min(nlx,nly) < 3)                                   /* Too few cells: use direct sum */
	{
		for (n1_sph = 1; n1_sph < num_sph; n1_sph++)
		for (n2_sph = 0; n2_sph < n1_sph;  n2_sph++)
		{

			x12 = n_image(particles[n1_sph].x0.x - particles[n2_sph].x0.x, demSize[0]);
			y12 = n_image(particles[n1_sph].x0.y - particles[n2_sph].x0.y, demSize[1]);
			z12 = n_image(particles[n1_sph].x0.z - particles[n2_sph].x0.z, demSize[2]);

			r12 = x12*x12 + y12*y12 + z12*z12;
			sig = particles[n1_sph].r + particles[n2_sph].r;
			cut = sig + 2.0;


			if (r12 <= cut*cut)
			{
				particles[n1_sph].list[0]++;
				particles[n2_sph].list[0]++;
				particles[n1_sph].xix[0]+= 1.0 ;
				particles[n1_sph].xiy[0]+= 1.0;
				particles[n1_sph].xiz[0]+= 1.0;
				particles[n2_sph].xix[0]+= 1.0;
				particles[n2_sph].xiy[0]+= 1.0;
				particles[n2_sph].xiz[0]+= 1.0;

				particles[n1_sph].list[particles[n1_sph].list[0]] = n2_sph;
				particles[n2_sph].list[particles[n2_sph].list[0]] = n1_sph;

				particles[n1_sph].xix[(int) particles[n1_sph].list[0]] =  particles[n2_sph].xix[(int) particles[n2_sph].list[0]];
				particles[n2_sph].xix[(int) particles[n2_sph].list[0]] =  particles[n1_sph].xix[(int) particles[n1_sph].list[0]];

				particles[n1_sph].xiy[(int) particles[n1_sph].list[0]] =  particles[n2_sph].xiy[(int) particles[n2_sph].list[0]];
				particles[n2_sph].xiy[(int) particles[n2_sph].list[0]] =  particles[n1_sph].xiy[(int) particles[n1_sph].list[0]];

				particles[n1_sph].xiz[(int) particles[n1_sph].list[0]] =  particles[n2_sph].xiz[(int) particles[n2_sph].list[0]];
				particles[n2_sph].xiz[(int) particles[n2_sph].list[0]] =  particles[n1_sph].xiz[(int) particles[n1_sph].list[0]];
			}
		}
	}
	else                                                               /* Use cell structure */
	{


		for (lx1 = 0; lx1 < nlx; lx1++)
		for (ly1 = 0; ly1 < nly; ly1++)
		for (lz1 = 0; lz1 < nlz; lz1++)
		{
			l1  = lx1*nly*nlz + ly1*nlz + lz1;
			cell_list[l1][0] = 0;
		}

		for (n1_sph = 0; n1_sph < num_sph; n1_sph++)
		{
			lx1 = (int) (box(particles[n1_sph].x0.x/demSize[0], 1)*nlx);
			ly1 = (int) (box(particles[n1_sph].x0.y/demSize[1], 1)*nly);
			lz1 = (int) (box(particles[n1_sph].x0.z/demSize[2], 1)*nlz);
			l1  = lx1*nly*nlz + ly1*nlz + lz1;

			cell_list[l1][0]++;

			cell_list[l1][cell_list[l1][0]] = n1_sph;

			particles[n1_sph].list[0] = 0;

			particles[n1_sph].xix[0] = 0.0;
	  	particles[n1_sph].xiy[0] = 0.0;
			particles[n1_sph].xiz[0] = 0.0;
		}

		for (lx1 = 0; lx1 < nlx; lx1++)
		for (ly1 = 0; ly1 < nly; ly1++)
		for (lz1 = 0; lz1 < nlz; lz1++)
		for (lx2 = lx1-1; lx2 <= lx1+1; lx2++)
		for (ly2 = ly1-1; ly2 <= ly1+1; ly2++)
		for (lz2 = lz1-1; lz2 <= lz1+1; lz2++)
		{
			l1  = lx1*nly*nlz + ly1*nlz + lz1;
			l2  = ((lx2+nlx)%nlx)*nly*nlz + ((ly2+nly)%nly)*nlz + (lz2+nlz)%nlz;
			for (n1 = 1; n1 <= cell_list[l1][0]; n1++)
			for (n2 = 1; n2 <= cell_list[l2][0]; n2++)
			{
				n1_sph  = cell_list[l1][n1];
				n2_sph  = cell_list[l2][n2];
				if (n2_sph < n1_sph)
				{
					x12 = n_image(particles[n1_sph].x0.x - particles[n2_sph].x0.x, demSize[0]);
					y12 = n_image(particles[n1_sph].x0.y - particles[n2_sph].x0.y, demSize[1]);
					z12 = n_image(particles[n1_sph].x0.z - particles[n2_sph].x0.z, demSize[2]);
					r12 = x12*x12 + y12*y12 + z12*z12;

					sig = particles[n1_sph].r + particles[n2_sph].r;
					cut = sig + 2.0;
					if (r12 <= cut*cut)
					{

						particles[n1_sph].list[0]++;
						particles[n2_sph].list[0]++;
						particles[n1_sph].xix[0]+= 1.0;
						particles[n1_sph].xiy[0]+= 1.0;
						particles[n1_sph].xiz[0]+= 1.0;
						particles[n2_sph].xix[0]+= 1.0;
						particles[n2_sph].xiy[0]+= 1.0;
						particles[n2_sph].xiz[0]+= 1.0;


						particles[n1_sph].list[particles[n1_sph].list[0]] = n2_sph;
						particles[n2_sph].list[particles[n2_sph].list[0]] = n1_sph;

					}
				}
			}
		}
	}




	for(n_sph = 0; n_sph < num_sph ; n_sph++)
	{
		for(n=1; n <= particles[n_sph].list[0];n++)
		{
			k = particles[n_sph].list[n];

			for (l=1; l <= particles[n_sph]._list[0]; l++)
			{
				j = particles[n_sph]._list[l];

				if (j == k)
				{
					particles[n_sph].xix[n] = particles[n_sph]._xix[l];
					particles[n_sph].xiy[n] = particles[n_sph]._xiy[l];
					particles[n_sph].xiz[n] = particles[n_sph]._xiz[l];
				}
			}
		}
	}
}

// evaluate all forces (including contact and hydrodynamics)
void DEM::evaluateForces(int timeStep){

	totalStress.reset();
	contactStress.reset();
	lubStress.reset();
	electroStress.reset();
	hydroStress.reset();
	//wallForce = 0.0;

		for (int n = 0; n < particles.size(); n++){

			particles[n].FTotal.reset();
			particles[n].MTotal.reset();
			particles[n].FContact.reset();
			particles[n].MContact.reset();
			particles[n].FLub.reset();
			particles[n].MLub.reset();
			particles[n].SLub.reset();
			particles[n].SContact.reset();
			particles[n].SRep.reset();
			particles[n].STotal.reset();
			particles[n].FRep.reset();
      particles[n].FElas.reset();
      particles[n].FVisc.reset();
      particles[n].SElas.reset();
      particles[n].SVisc.reset();
      particles[n].FVisc_h.reset();
      particles[n].SVisc_h.reset();
      particles[n].coord = 0.0;
      particles[n].lub_coord = 0.0;

		}

		for (int o = 0; o < objects.size(); o++) {

			objects[o].FTotal.reset();
			objects[o].FContact.reset();
			objects[o].FLub.reset();
      objects[o].FSpring.reset();
      objects[o].MTotal.reset();
      objects[o].MContact.reset();
      objects[o].MLub.reset();
      objects[o].FElas.reset();
      objects[o].FVisc.reset();
      objects[o].FVisc_h.reset();

		}

    for (int w = 0; w < walls.size(); w++) {
 			walls[w].FTotal.reset();
 			walls[w].FContact.reset();
 			walls[w].FLub.reset();
			walls[w].FRep.reset();
		 }

     for (int b = 0; b < bodies.size(); b++) {
 			bodies[b].FTotal.reset();
      bodies[b].FSpring.reset();
 		}


		if (demTimeStep >= objWaitTime){

      if (bodies.size()){
        bodyObjectInteraction();
      }

      calculateParticleLubricationForces();
      calculateWallLubricationForces();
      calculateObjectParticleContactForces();
      calculateObjectParticleLubricationForces();
  } else{

    calculateFictionalWallContactForces();

  }


		calculateParticleContactForces();


		calculateParticleElectroForces();


		calculateWallContactForces();





		for (int n = 0; n < particles.size(); n++) {
			lubStress += particles[n].SLub;
			hydroStress += particles[n].SHydro;
      particles[n].STotal += particles[n].SHydro + particles[n].SLub + particles[n].SContact + particles[n].SRep ;
      particles[n].SVisc_h += particles[n].SHydro + particles[n].SLub;
      particles[n].SVisc += particles[n].SVisc_h;
      particles[n].pressure = -1.0/3.0 * (particles[n].STotal.m00 + particles[n].STotal.m11 + particles[n].STotal.m22);
      particles[n].normal1 = -particles[n].STotal.m00 + particles[n].STotal.m11;
      particles[n].normal2 = -particles[n].STotal.m11 + particles[n].STotal.m22;
      particles[n].mu = (-particles[n].STotal.m12) / particles[n].pressure;
      particles[n].ratio = (-particles[n].STotal.m11)/ (-particles[n].STotal.m12);
      particles[n].trueMod = particles[n].STotal.m11 / particles[n].trueStrain;
      particles[n].engMod = particles[n].STotal.m11 / particles[n].engStrain;
      particles[n].FTotal += particles[n].FHydro + particles[n].FContact + particles[n].FLub + particles[n].FGrav + particles[n].FRep;
      particles[n].MTotal += particles[n].MHydro + particles[n].MContact + particles[n].MLub;
      particles[n].FVisc_h += particles[n].FHydro + particles[n].FLub;
      particles[n].overlap = particles[n].FElas / sphereMat.linearStiff;
		}



		lubStress /= vol;
		contactStress /= vol;
		electroStress /= vol;
		hydroStress /= vol;

		totalStress += contactStress + lubStress + electroStress + hydroStress;

		for (int w = 0; w < walls.size(); w++) {
			walls[w].FTotal += walls[w].FHydro + walls[w].FContact + walls[w].FLub + walls[w].FRep;
    	}

    if (demTimeStep <= objWaitTime){

      for (int n = 0; n < particles.size(); n++){
        	particles[n].FGrav = tVect(0.0,0.0,0.0);
      }

    }


		for (int o = 0; o < objects.size(); o++) {
			objects[o].FTotal += objects[o].FHydro + objects[o].FContact + objects[o].FLub + objects[o].FGrav + objects[o].FSpring;
      objects[o].MTotal += objects[o].MHydro + objects[o].MContact + objects[o].MLub;
      objects[o].FVisc_h += objects[o].FHydro + objects[o].FLub;
		}

    for (int o = 0; o < bodies.size(); o++) {
      bodies[o].FTotal += bodies[o].FGrav + bodies[o].FSpring;
    }

}

void DEM::bodyObjectInteraction(){

for (int b = 0; b < bodies.size(); b++){
  for (int o = 0; o < objects.size(); o++){
    double fs;
    double fdo;
    double fdb;
    bodies[b].L = bodies[b].x0.y - objects[o].x0.y;
  //  double avgM =  0.5*(objects[o].m + bodies[b].m);
    double avgM = objects[o].m ;
    tVect r = bodies[b].x0 - objects[o].x0;
    tVect n = r / r.norm();
    tVect ur = objects[o].x1 - bodies[b].x1;
    double gn = ur.x*n.x + ur.z*n.z + ur.y*n.y;
  //  double bodyDamp = 0.0;
    double del = bodies[b].L0- bodies[b].L;

    double bodyDamp = 2.0*sqrt(bodies[b].m*bodies[b].k) * springDamp;
    double objDamp  = 2.0*sqrt(objects[o].m*bodies[b].k) * springDamp;

    double sigma = bodies[b].r ;

    if (nonLinType == 0){
      // spring
      fs = bodies[b].k * del*n.y;
      fdo = 0.0;
      fdb = 0.0;

    }
    else if (nonLinType == 1){
      // spring-dashpot
      fs = (bodies[b].k * del )*n.y;
      fdo = -objects[o].x1.y * objDamp;
      fdb = -bodies[b].x1.y * bodyDamp;

    }
    else if (nonLinType == 2){
      //spring-dashpot-cubic
      fs = (bodies[b].k * del + k13*bodies[b].k*(del*del*del))*n.y;
      fdo = -objects[o].x1.y * objDamp;
      fdb = -bodies[b].x1.y * bodyDamp;

    }
    else if (nonLinType == 3){
      //spring-dashpot-repLJ
      fs = (bodies[b].k * del  + (24.0*epsilon/r.norm())*2.0*(pow((sigma/r.norm()), 12.0)))*n.y;
      fdo = -objects[o].x1.y * objDamp;
      fdb = -bodies[b].x1.y * bodyDamp;

    }
    else if (nonLinType == 4){
      //spring-dashpot-WCA
      double wcacut = pow(2.0,1.0/6.0)*sigma;
      if (r.norm() <= wcacut ){
      fs = (bodies[b].k * del  + (24.0*epsilon/r.norm()) * ( ( 2.0* pow( sigma/r.norm() , 12.0) )  - ( pow( sigma/r.norm() , 6.0) )  ))*n.y;
      fdo = -objects[o].x1.y * objDamp;
      fdb = -bodies[b].x1.y * bodyDamp;
      }
      else{
      fs = (bodies[b].k * del )*n.y;
      fdo = -objects[o].x1.y * objDamp;
      fdb = -bodies[b].x1.y * bodyDamp;
      }
    }

    objects[o].FSpring.y -= fs;
    bodies[b].FSpring.y += fs;
    objects[o].FSpring.y += fdo;
    bodies[b].FSpring.y += fdb;
    }
  }
}


tVect DEM::calculateEdgeDistance(tVect rb, double w2, double h2, double d2){
  double s;
  tVect rbs = tVect(abs(rb.x)/w2, abs(rb.y)/h2, abs(rb.z)/d2);

  if (rbs.x == rbs.y == rbs.z){
    //corners
    double diag = sqrt(w2*w2 + h2*h2 + d2*d2);
    s = diag / rb.norm();
  }
  else if (rbs.x == rbs.y && rbs.x != rbs.z){
    // xy-intersect
    double diag = sqrt(w2*w2 + h2*h2);
    s = diag / sqrt(rb.x*rb.x + rb.y*rb.y);
  }
  else if (rbs.z == rbs.y && rbs.z != rbs.x){
    // yz-intersect
    double diag = sqrt(d2*d2 + h2*h2);
    s = diag / sqrt(rb.z*rb.z + rb.y*rb.y);
  }
  else if (rbs.z == rbs.x && rbs.z != rbs.y){
    // xz-intersect
    double diag = sqrt(d2*d2 + w2*w2);
    s = diag / sqrt(rb.z*rb.z + rb.x*rb.x);
  }
  else if (max(rbs.y,max(rbs.x,rbs.z)) == rbs.x){
    // x-intersect
    s = w2 / abs(rb.x);
  }
  else if (max(rbs.y,max(rbs.x,rbs.z)) == rbs.y){
    // y-intersect
    s = h2 /abs(rb.y);
  }
  else if (max(rbs.y,max(rbs.x,rbs.z)) == rbs.z){
    // z-intersect
    s = d2 / abs(rb.z);
  }

  return rb * s;
}

void DEM::calculateObjectParticleContactForces() {

	// object-particle contact
  tVect cpoint, cpoint_w;
	double o1wx,o1wy,o1wz, v1wx, v1wy, v1wz, xi1w, fnxw, fnyw,fnzw, ktw, dampCoeffTanw, ftxw, ftyw, ftzw;
	double fnnw, fttw, friw, t1wx, t1wy, t1wz, dfxw, dfyw, dfzw, dtxw, dtyw, dtzw;
	double n1wx, n1wy, n1wz, gn1w, muw, radw, delw, r1w, dis;
  double felx, fely, felz, fvisx, fvisy, fvisz;

	for (int n1 = 0; n1 < particles.size(); n1++){

		for (int n2_w = 0; n2_w < objects.size(); n2_w++){

      double radw = 0.5*(particles[n1].r + objects[n2_w].r) ;
      double x1w =  particles[n1].x0.x - objects[n2_w].x0.x ;
      double y1w =  particles[n1].x0.y - objects[n2_w].x0.y ;
      double z1w =  particles[n1].x0.z - objects[n2_w].x0.z ;
      double r1w = sqrt(x1w*x1w + y1w*y1w + z1w*z1w);

      double n1wx = x1w / r1w;
      double n1wy = y1w / r1w;
      double n1wz = z1w / r1w;


    double urxw = particles[n1].x1.x - objects[n2_w].x1.x;
    double uryw = particles[n1].x1.y - objects[n2_w].x1.y;
    double urzw = particles[n1].x1.z - objects[n2_w].x1.z;

    double delw = (2*radw)- r1w;
    //  double zeta = sqrt(sphereMat.linearStiff*((particles[n1].m + objects[n2_w].m)*0.5));

			 if (delw > 0.0){

			 	gn1w = urxw*n1wx + uryw*n1wy + urzw*n1wz;

				muw = sphereMat.frictionCoefPart;

			  if (muw != 0.0 ){

				o1wx = particles[n1].r * particles[n1].w1.x  + objects[n2_w].r * objects[n2_w].w1.x;
				o1wy = particles[n1].r * particles[n1].w1.y  + objects[n2_w].r  * objects[n2_w].w1.y;
				o1wz = particles[n1].r * particles[n1].w1.z  + objects[n2_w].r * objects[n2_w].w1.z;

				v1wx = urxw + n1wy*o1wz - n1wz*o1wy - gn1w*n1wx;
				v1wy = uryw + n1wz*o1wx - n1wx*o1wz - gn1w*n1wy;
				v1wz = urzw + n1wx*o1wy - n1wy*o1wx - gn1w*n1wz;

				xi1w = particles[n1].xio.x*n1wx + particles[n1].xio.y*n1wy + particles[n1].xio.z*n1wz ;

				particles[n1].xio.x -= xi1w*n1wx;
				particles[n1].xio.y -= xi1w*n1wy;
				particles[n1].xio.z -= xi1w*n1wz;


				fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
				fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
				fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;

        felx = sphereMat.linearStiff*delw*n1wx;
        fely = sphereMat.linearStiff*delw*n1wy;
        felz = sphereMat.linearStiff*delw*n1wz;

        fvisx = - sphereMat.dampCoeff * gn1w * n1wx;
        fvisy = - sphereMat.dampCoeff * gn1w * n1wy;
        fvisz = - sphereMat.dampCoeff * gn1w * n1wz;
				//tangential stiffness
				ktw = 0.2*sphereMat.linearStiff;


				ftxw = - ktw * particles[n1].xio.x ;
				ftyw = - ktw * particles[n1].xio.y ;
				ftzw = - ktw * particles[n1].xio.z ;

				fnnw = sqrt(fnxw*fnxw + fnyw*fnyw + fnzw*fnzw);
				fttw = sqrt(ftxw*ftxw + ftyw*ftyw + ftzw*ftzw);

				friw = muw*fnnw;

				t1wx = ftxw/fttw;
				t1wy = ftyw/fttw;
				t1wz = ftzw/fttw;


				if (friw < fttw){

					ftxw = friw*t1wx;
					ftyw = friw*t1wy;
					ftzw = friw*t1wz;

          double fne = sqrt(felx*felx + fely*fely + felz*felz);
          double fnv = sqrt(fvisx*fvisx + fvisy*fvisy + fvisz*fvisz);

          felx += muw*fne*t1wx;
          fely += muw*fne*t1wy;
          fely += muw*fne*t1wz;

          fvisx += muw*fnv*t1wx;
          fvisx += muw*fnv*t1wy;
          fvisx += muw*fnv*t1wz;


					particles[n1].xio.x = - (1/ktw) * (ftxw + dampCoeffTanw * v1wx) ;
					particles[n1].xio.y = - (1/ktw) * (ftyw + dampCoeffTanw * v1wy) ;
					particles[n1].xio.z = - (1/ktw) * (ftzw + dampCoeffTanw * v1wz) ;

				}
        else{
          felx +=   - ktw * particles[n1].xio.x;
          fely +=   - ktw * particles[n1].xio.y;
          felz +=   - ktw * particles[n1].xio.z;

          fvisx += - dampCoeffTanw * v1wx;
          fvisy += - dampCoeffTanw * v1wy;
          fvisz += - dampCoeffTanw * v1wz;

          particles[n1].xio.x += v1wx*deltat;
          particles[n1].xio.y += v1wy*deltat;
          particles[n1].xio.z += v1wz*deltat;
        }

			}
				else {
					 fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
					 fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
					 fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;

           felx = sphereMat.linearStiff*delw*n1wx;
           fely = sphereMat.linearStiff*delw*n1wy;
           felz = sphereMat.linearStiff*delw*n1wz;

           fvisx = - sphereMat.dampCoeff * gn1w * n1wx;
           fvisy = - sphereMat.dampCoeff * gn1w * n1wy;
           fvisz = - sphereMat.dampCoeff * gn1w * n1wz;

					 ftxw = 0.0;
					 ftyw = 0.0;
					 ftzw = 0.0;

				}
			}
			else {

				fnxw = fnyw = fnzw = ftxw = ftyw = ftzw = 0.0;
        felx = fely = felz = fvisx = fvisy = fvisz =0.0;

				particles[n1].xio.x = 0.0;
				particles[n1].xio.y = 0.0;
				particles[n1].xio.z = 0.0;

			}

			double dfxw = fnxw + ftxw;
			double dfyw = fnyw + ftyw;
			double dfzw = fnzw + ftzw;

			double dtxw = ftyw*n1wz - ftzw*n1wy;
			double dtyw = ftzw*n1wx - ftxw*n1wz;
			double dtzw = ftxw*n1wy - ftyw*n1wx;


      tVect totForce = tVect(dfxw,dfyw,dfzw);

			particles[n1].FContact.x += dfxw;
			particles[n1].FContact.y += dfyw;
			particles[n1].FContact.z += dfzw;

			objects[n2_w].FContact.x -= dfxw;
			objects[n2_w].FContact.y -= dfyw;
			objects[n2_w].FContact.z -= dfzw;

			particles[n1].MContact.x += particles[n1].r * dtxw ;
			particles[n1].MContact.y += particles[n1].r * dtyw ;
			particles[n1].MContact.z += particles[n1].r * dtzw ;


      objects[n2_w].MContact.x += objects[n2_w].r * dtxw ;
      objects[n2_w].MContact.y += objects[n2_w].r * dtyw ;
      objects[n2_w].MContact.z += objects[n2_w].r * dtzw ;

      particles[n1].FElas = particles[n1].FElas + tVect(felx,fely,felz);
      particles[n1].FVisc = particles[n1].FVisc + tVect(fvisx,fvisy,fvisz);
      objects[n2_w].FElas = objects[n2_w].FElas - tVect(felx,fely,felz);
      objects[n2_w].FVisc = objects[n2_w].FVisc - tVect(fvisx,fvisy,fvisz);


			particles[n1].SContact.m00 += -x1w  * dfxw;
			particles[n1].SContact.m01 += -x1w  * dfyw;
			particles[n1].SContact.m02 += -x1w  * dfzw;
			particles[n1].SContact.m10 += -y1w  * dfxw;
			particles[n1].SContact.m11 += -y1w  * dfyw;
			particles[n1].SContact.m12 += -y1w  * dfzw;
			particles[n1].SContact.m20 += -z1w  * dfxw;
			particles[n1].SContact.m21 += -z1w  * dfyw;
			particles[n1].SContact.m22 += -z1w  * dfzw;

      particles[n1].SElas.m00 += -x1w  * felx;
      particles[n1].SElas.m01 += -x1w  * fely;
      particles[n1].SElas.m02 += -x1w  * felz;
      particles[n1].SElas.m10 += -y1w  * felx;
      particles[n1].SElas.m11 += -y1w  * fely;
      particles[n1].SElas.m12 += -y1w  * felz;
      particles[n1].SElas.m20 += -z1w  * felx;
      particles[n1].SElas.m21 += -z1w  * fely;
      particles[n1].SElas.m22 += -z1w  * felz;

      particles[n1].SVisc.m00 += -x1w  * fvisx;
      particles[n1].SVisc.m01 += -x1w  * fvisy;
      particles[n1].SVisc.m02 += -x1w  * fvisz;
      particles[n1].SVisc.m10 += -y1w  * fvisx;
      particles[n1].SVisc.m11 += -y1w  * fvisy;
      particles[n1].SVisc.m12 += -y1w  * fvisz;
      particles[n1].SVisc.m20 += -z1w  * fvisx;
      particles[n1].SVisc.m21 += -z1w  * fvisy;
      particles[n1].SVisc.m22 += -z1w  * fvisz;

			contactStress.m00 += -x1w  * dfxw;
			contactStress.m01 += -x1w  * dfyw;
			contactStress.m02 += -x1w  * dfzw;
			contactStress.m10 += -y1w  * dfxw;
			contactStress.m11 += -y1w  * dfyw;
			contactStress.m12 += -y1w  * dfzw;
			contactStress.m20 += -z1w  * dfxw;
			contactStress.m21 += -z1w  * dfyw;
			contactStress.m22 += -z1w  * dfzw;
		}
	}
}

void DEM::calculateObjectParticleLubricationForces() {
// Also from Susp3D, based on Nguyen and Ladd, PRE (2002)

	double        Xa, Ya, Yb1, Yb2, Xg1, Xg2, Yg1, Yg2;
	double        Yc11, Yc12, Yc21, Yc22, Yh11, Yh12, Yh21, Yh22;
	double        Fx, Fy, Fz, T1x, T1y, T1z, T2x, T2y, T2z;
	double        S1xx, S1yy, S1zz, S1yz, S1zx, S1xy;
	double        S2xx, S2yy, S2zz, S2yz, S2zx, S2xy;
	double        nu, rad, cut, dr, h, h_c;
	double        umx, umy, umz, w1x, w1y, w1z, w2x, w2y, w2z;
	double        umxrx, umxry, umxrz, w1xrx, w1xry, w1xrz, w2xrx, w2xry, w2xrz;
	double        um_r, w1_r, w2_r, beta;

	double lub_N = sphereMat.lubNormal;
	double lub_T = sphereMat.lubTangential;
	double lub_R = sphereMat.lubRotation;


	for (int n1 =0; n1 < particles.size(); n1++){

		for (int n2 = 0; n2 < objects.size(); n2++){

				cut = sphereMat.lubCutoff;
				rad = 0.5*(particles[n1].r + objects[n2].r);
				double x12 = particles[n1].x0.x - objects[n2].x0.x ;
				double y12 = particles[n1].x0.y - objects[n2].x0.y ;
				double z12 = particles[n1].x0.z - objects[n2].x0.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

				umx = particles[n1].x1.x - objects[n2].x1.x;
				umy = particles[n1].x1.y - objects[n2].x1.y;
				umz = particles[n1].x1.z - objects[n2].x1.z;

				w1x   = particles[n1].w1.x;
				w1y   = particles[n1].w1.y;
				w1z   = particles[n1].w1.z;
				w2x   = objects[n2].w1.x;
				w2y   = objects[n2].w1.y;
				w2z   = objects[n2].w1.z;

				double dr = r12 - 2.0*rad;
				double h = max(dr,cut);

				if ( dr < max(cut, max(lub_N, max(lub_T,lub_R))))
				{

					double nu = sphereMat.kinVisc;
					double alpha = M_PI*nu*rad*deltat;
					double beta1 = objects[n2].r / particles[n1].r;
					double beta2 = beta1*beta1;
					double beta3 = beta1*beta2;
					double beta4 = pow((1+beta1), 4.0);
					double beta5 = pow((1+beta1), 5.0);

					if (dr < max(cut, lub_N))
					{
						Xa   =  4.0*beta2/beta4;                         /* Set lubrication coefficients */
						Xg1  = 12.0*beta2/beta5;
						Xg2  =-12.0*beta3/beta5;
						h_c   = max(dr, lub_N);
						Xa   *= rad/h - rad/h_c;
						Xg1  *= rad/h - rad/h_c;
						Xg2  *= rad/h - rad/h_c;
						Xa   *= 6.0*alpha;
						Xg1  *= 4.0*alpha*rad;
						Xg2  *= 4.0*alpha*rad;
					}
					else
						Xa = Xg1 = Xg2 = 0.0;
					if (dr < max(cut, lub_T))
					{
						Ya    = 8.0*beta1*(2.0+beta1+2.0*beta2)/(15.0*beta4);
						Yb1   =-4.0*beta1*(4.0+beta1)/(5.0*beta4);
						Yb2   = 4.0*beta1*(4.0*beta2+beta1)/(5.0*beta4);
						Yg1   = 2.0*beta1*(4.0-beta1+7.0*beta2)/(5.0*beta5);
						Yg2   =-2.0*beta2*(4.0*beta2-beta1+7.0)/(5.0*beta5);
						h_c   = max(dr, lub_T);
						beta  = alpha*log(h_c/h);
						Ya   *= 6.0*beta;
						Yb1  *= 4.0*beta*rad;
						Yb2  *= 4.0*beta*rad;
						Yg1  *= 4.0*beta*rad;
						Yg2  *= 4.0*beta*rad;
					}
					else
						Ya = Yb1 = Yb2 = Yg1 = Yg2 = 0.0;
					if (dr < max(cut, lub_R))
					{
						Yc11  = 16.0*beta1/(5.0*beta4);
						Yc12  =  4.0*beta2/(5.0*beta4);
						Yc21  =  4.0*beta2/(5.0*beta4);
						Yc22  = 16.0*beta3/(5.0*beta4);
						Yh11  =  4.0*beta1*(2.0-beta1)/(5.0*beta5);
						Yh12  =  2.0*beta2*(1.0+7.0*beta1)/(5.0*beta5);
						Yh21  =  2.0*beta2*(beta1+7.0)/(5.0*beta5);
						Yh22  =  4.0*beta2*(2.0*beta2-beta1)/(5.0*beta5);
						h_c   = max(dr, lub_R);
						beta  = alpha*log(h_c/h);
						Yc11 *= 8.0*beta*rad*rad;
						Yc12 *= 8.0*beta*rad*rad;
						Yc21 *= 8.0*beta*rad*rad;
						Yc22 *= 8.0*beta*rad*rad;
						Yh11 *= 8.0*beta*rad*rad;
						Yh21 *= 8.0*beta*rad*rad;
						Yh12 *= 8.0*beta*rad*rad;
						Yh22 *= 8.0*beta*rad*rad;
					}
					else
						Yc11 = Yc12 = Yc21 = Yc22 = Yh11 = Yh12 = Yh21 = Yh22 = 0.0;

					x12  /= r12;
					y12  /= r12;
					z12  /= r12;

					um_r  = umx*x12 + umy*y12 + umz*z12;
					w1_r  = w1x*x12 + w1y*y12 + w1z*z12;
					w2_r  = w2x*x12 + w2y*y12 + w2z*z12;
					umxrx = umy*z12 - umz*y12;
					umxry = umz*x12 - umx*z12;
					umxrz = umx*y12 - umy*x12;
					w1xrx = w1y*z12 - w1z*y12;
					w1xry = w1z*x12 - w1x*z12;
					w1xrz = w1x*y12 - w1y*x12;
					w2xrx = w2y*z12 - w2z*y12;
					w2xry = w2z*x12 - w2x*z12;
					w2xrz = w2x*y12 - w2y*x12;

					Fx    = -Xa*um_r*x12 - Ya*(umx-um_r*x12);
					Fy    = -Xa*um_r*y12 - Ya*(umy-um_r*y12);
					Fz    = -Xa*um_r*z12 - Ya*(umz-um_r*z12);
					Fx   += -Yb1*w1xrx + Yb2*w2xrx;
					Fy   += -Yb1*w1xry + Yb2*w2xry;
					Fz   += -Yb1*w1xrz + Yb2*w2xrz;
					T1x   =  Yb1*umxrx;
					T1y   =  Yb1*umxry;
					T1z   =  Yb1*umxrz;
					T2x   = -Yb2*umxrx;
					T2y   = -Yb2*umxry;
					T2z   = -Yb2*umxrz;
					T1x  += -Yc11*(w1x-w1_r*x12) - Yc12*(w2x-w2_r*x12);
					T1y  += -Yc11*(w1y-w1_r*y12) - Yc12*(w2y-w2_r*y12);
					T1z  += -Yc11*(w1z-w1_r*z12) - Yc12*(w2z-w2_r*z12);
					T2x  += -Yc21*(w1x-w1_r*x12) - Yc22*(w2x-w2_r*x12);
					T2y  += -Yc21*(w1y-w1_r*y12) - Yc22*(w2y-w2_r*y12);
					T2z  += -Yc21*(w1z-w1_r*z12) - Yc22*(w2z-w2_r*z12);
					S1xx  =  Xg1*um_r*x12*x12 + Yg1*(x12*umx+umx*x12);
					S1yy  =  Xg1*um_r*y12*y12 + Yg1*(y12*umy+umy*y12);
					S1zz  =  Xg1*um_r*z12*z12 + Yg1*(z12*umz+umz*z12);
					S1yz  =  Xg1*um_r*y12*z12 + Yg1*(y12*umz+umy*z12);
					S1zx  =  Xg1*um_r*z12*x12 + Yg1*(z12*umx+umz*x12);
					S1xy  =  Xg1*um_r*x12*y12 + Yg1*(x12*umy+umx*y12);
					S2xx  = -Xg2*um_r*x12*x12 - Yg2*(x12*umx+umx*x12);
					S2yy  = -Xg2*um_r*y12*y12 - Yg2*(y12*umy+umy*y12);
					S2zz  = -Xg2*um_r*z12*z12 - Yg2*(z12*umz+umz*z12);
					S2yz  = -Xg2*um_r*y12*z12 - Yg2*(y12*umz+umy*z12);
					S2zx  = -Xg2*um_r*z12*x12 - Yg2*(z12*umx+umz*x12);
					S2xy  = -Xg2*um_r*x12*y12 - Yg2*(x12*umy+umx*y12);
					S1xx += -Yh11*(x12*w1xrx+w1xrx*x12) - Yh12*(x12*w2xrx+w2xrx*x12);
					S1yy += -Yh11*(y12*w1xry+w1xry*y12) - Yh12*(y12*w2xry+w2xry*y12);
					S1zz += -Yh11*(z12*w1xrz+w1xrz*z12) - Yh12*(z12*w2xrz+w2xrz*z12);
					S1yz += -Yh11*(y12*w1xrz+w1xry*z12) - Yh12*(y12*w2xrz+w2xry*z12);
					S1zx += -Yh11*(z12*w1xrx+w1xrz*x12) - Yh12*(z12*w2xrx+w2xrz*x12);
					S1xy += -Yh11*(x12*w1xry+w1xrx*y12) - Yh12*(x12*w2xry+w2xrx*y12);
					S2xx += -Yh21*(x12*w1xrx+w1xrx*x12) - Yh22*(x12*w2xrx+w2xrx*x12);
					S2yy += -Yh21*(y12*w1xry+w1xry*y12) - Yh22*(y12*w2xry+w2xry*y12);
					S2zz += -Yh21*(z12*w1xrz+w1xrz*z12) - Yh22*(z12*w2xrz+w2xrz*z12);
					S2yz += -Yh21*(y12*w1xrz+w1xry*z12) - Yh22*(y12*w2xrz+w2xry*z12);
					S2zx += -Yh21*(z12*w1xrx+w1xrz*x12) - Yh22*(z12*w2xrx+w2xrz*x12);
					S2xy += -Yh21*(x12*w1xry+w1xrx*y12) - Yh22*(x12*w2xry+w2xrx*y12);


					particles[n1].FLub.x += Fx;
					particles[n1].FLub.y += Fy;
					particles[n1].FLub.z += Fz;
					objects[n2].FLub.x -= Fx;
					objects[n2].FLub.y -= Fy;
					objects[n2].FLub.z -= Fz;

					particles[n1].MLub.x += T1x;
					particles[n1].MLub.y += T1y;
					particles[n1].MLub.z += T1z;
					objects[n2].MLub.x += T2x;
					objects[n2].MLub.y += T2y;
					objects[n2].MLub.z += T2z;

					particles[n1].SLub.m00 += S1xx;
					particles[n1].SLub.m01 += S1xy;
					particles[n1].SLub.m02 += S1zx;
					particles[n1].SLub.m10 += S1xy;
					particles[n1].SLub.m11 += S1yy;
					particles[n1].SLub.m12 += S1yz;
					particles[n1].SLub.m20 += S1zx;
					particles[n1].SLub.m21 += S1yz;
					particles[n1].SLub.m22 += S1zz;


			}
		}
	}

}


void DEM::calculateWallElectroForces(){

	for (int n1 =0 ; n1 < particles.size(); n1++){

		for (int n2 = 0; n2 < walls.size(); n2++){


				double rad = particles[n1].r;
				double x12 = (particles[n1].x0.x - walls[n2].p.x) * walls[n2].n.x * walls[n2].n.x;
				double y12 = (particles[n1].x0.y - walls[n2].p.y) * walls[n2].n.y * walls[n2].n.y;
				double z12 = (particles[n1].x0.z - walls[n2].p.z) * walls[n2].n.z * walls[n2].n.z;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);


				double	n12x = x12/r12;
				double	n12y = y12/r12;
				double	n12z = z12/r12;


				double fer = sphereMat.electroForce;
				double lamda = sphereMat.debyeLength;
				double del = rad - r12;

				double frx,fry,frz;

				if ( del <= 0.0){
					frx= fer*exp(del/lamda)*n12x;
					fry= fer*exp(del/lamda)*n12y;
					frz= fer*exp(del/lamda)*n12z;

				}
				else {
					frx= fer*n12x;
					fry= fer*n12y;
					frz= fer*n12z;
				}

					particles[n1].FRep.x += frx;
					particles[n1].FRep.y += fry;
					particles[n1].FRep.z += frz;

					walls[n2].FRep.x -= frx;
					walls[n2].FRep.y -= fry;
					walls[n2].FRep.z -= frz;

					electroStress.m00 += -x12  * frx;
					electroStress.m01 += -x12  * fry;
					electroStress.m02 += -x12  * frz;
					electroStress.m10 += -y12  * frx;
					electroStress.m11 += -y12  * fry;
					electroStress.m12 += -y12  * frz;
					electroStress.m20 += -z12  * frx;
					electroStress.m21 += -z12  * fry;
					electroStress.m22 += -z12  * frz;


		}
	}

}

void DEM::calculateParticleElectroForces(){

	for (int n1 =0 ; n1 < particles.size(); n1++){

		for (int n = 1; n <= particles[n1].list[0] ; n++){

			int n2 = particles[n1].list[n];


			if (n2 < n1){

				double rad = 0.5*(particles[n1].r + particles[n2].r);
				double x12 = particles[n1].x0.x - particles[n2].x0.x ;
				double y12 = particles[n1].x0.y - particles[n2].x0.y ;
				double z12 = particles[n1].x0.z - particles[n2].x0.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);


				double	n12x = x12/r12;
				double	n12y = y12/r12;
				double	n12z = z12/r12;


				double fer = sphereMat.electroForce;
				double lamda = sphereMat.debyeLength;
				double del = 2*rad - r12;
				double frx,fry,frz;

				if ( del <= 0.0){
					frx= fer*exp(del/lamda)*n12x;
					fry= fer*exp(del/lamda)*n12y;
					frz= fer*exp(del/lamda)*n12z;

				}
				else {
					frx=fer*n12x;
					fry=fer*n12y;
					frz=fer*n12z;
				}

					particles[n1].FRep.x += frx;
					particles[n1].FRep.y += fry;
					particles[n1].FRep.z += frz;
					particles[n2].FRep.x -= frx;
					particles[n2].FRep.y -= fry;
					particles[n2].FRep.z -= frz;

					particles[n1].SRep.m00 += -x12  * frx;
					particles[n1].SRep.m01 += -x12  * fry;
					particles[n1].SRep.m02 += -x12  * frz;
					particles[n1].SRep.m10 += -y12  * frx;
					particles[n1].SRep.m11 += -y12  * fry;
					particles[n1].SRep.m12 += -y12  * frz;
					particles[n1].SRep.m20 += -z12  * frx;
					particles[n1].SRep.m21 += -z12  * fry;
					particles[n1].SRep.m22 += -z12  * frz;

					particles[n2].SRep.m00 += -x12  * frx;
					particles[n2].SRep.m01 += -x12  * fry;
					particles[n2].SRep.m02 += -x12  * frz;
					particles[n2].SRep.m10 += -y12  * frx;
					particles[n2].SRep.m11 += -y12  * fry;
					particles[n2].SRep.m12 += -y12  * frz;
					particles[n2].SRep.m20 += -z12  * frx;
					particles[n2].SRep.m21 += -z12  * fry;
					particles[n2].SRep.m22 += -z12  * frz;

					electroStress.m00 += -x12  * frx;
					electroStress.m01 += -x12  * fry;
					electroStress.m02 += -x12  * frz;
					electroStress.m10 += -y12  * frx;
					electroStress.m11 += -y12  * fry;
					electroStress.m12 += -y12  * frz;
					electroStress.m20 += -z12  * frx;
					electroStress.m21 += -z12  * fry;
					electroStress.m22 += -z12  * frz;

			}
		}
	}
}

// Calculate contact forces on particles
void DEM::calculateParticleContactForces() {
	// particle-particle contacts, from 3D-DEM by Satoshi Takada
	double o12x,o12y,o12z, v12x, v12y, v12z, xi12, fnx, fny,fnz, kt, dampCoeffTan, ftx, fty, ftz;
	double fnn, ftt, fri, t12x, t12y, t12z, dfx, dfy, dfz, dtx, dty, dtz;
	double x1, x2,y1,y2,z1,z2;
  double felx, fely, felz, fvisx, fvisy, fvisz;

	numContact = 0;

	numTriangle = 0;

	avgForce = 0.0;

	for (int n1 =0 ; n1 < particles.size(); n1++) {

		for (int n = 1; n <= particles[n1].list[0] ; n++){

			int n2 = particles[n1].list[n];

			if (n2 < n1){

				double rad = 0.5*(particles[n1].r + particles[n2].r);
				double x12 = particles[n1].x0.x - particles[n2].x0.x ;
				double y12 = particles[n1].x0.y - particles[n2].x0.y ;
				double z12 = particles[n1].x0.z - particles[n2].x0.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

				double urx = particles[n1].x1.x - particles[n2].x1.x;
				double ury = particles[n1].x1.y - particles[n2].x1.y;
				double urz = particles[n1].x1.z - particles[n2].x1.z;

				double	n12x = x12/r12;
				double	n12y = y12/r12;
				double	n12z = z12/r12;


				double gn12 = urx*n12x + ury*n12y + urz*n12z;

				double mu = sphereMat.frictionCoefPart;

				double del = 2*rad - r12;
				double dr = -del;


				if (del > 0.0) {

				if (mu != 0.0 ){


					o12x = particles[n1].r * particles[n1].w1.x + particles[n2].r * particles[n2].w1.x ;
					o12y = particles[n1].r * particles[n1].w1.y + particles[n2].r * particles[n2].w1.y ;
					o12z = particles[n1].r * particles[n1].w1.z + particles[n2].r * particles[n2].w1.z ;

					v12x = urx + n12y*o12z - n12z*o12y - gn12*n12x;
					v12y = ury + n12z*o12x - n12x*o12z - gn12*n12y;
				  v12z = urz + n12x*o12y - n12y*o12x - gn12*n12z;

				  xi12 = particles[n1].xix[n]*n12x + particles[n1].xiy[n]*n12y + particles[n1].xiz[n]*n12z ;

					particles[n1].xix[n] -= xi12*n12x;
					particles[n1].xiy[n] -= xi12*n12y;
					particles[n1].xiz[n] -= xi12*n12z;

				  fnx = sphereMat.linearStiff*del*n12x - sphereMat.dampCoeff * gn12 * n12x ;
					fny = sphereMat.linearStiff*del*n12y - sphereMat.dampCoeff * gn12 * n12y ;
					fnz = sphereMat.linearStiff*del*n12z - sphereMat.dampCoeff * gn12 * n12z ;


          felx = sphereMat.linearStiff*del*n12x;
          fely = sphereMat.linearStiff*del*n12y;
          felz = sphereMat.linearStiff*del*n12z;

          fvisx = - sphereMat.dampCoeff  * gn12 * n12x;
          fvisy = - sphereMat.dampCoeff  * gn12 * n12y;
          fvisz = - sphereMat.dampCoeff  * gn12 * n12z;

					//tangential stiffness
				  kt = 0.2*sphereMat.linearStiff;


					if (sphereMat.kinVisc == 0.0){
					dampCoeffTan = sphereMat.dampCoeff;
					}
					else {
					dampCoeffTan = 0.0; // no tangential dissipation as in LF-DEM by Seto and Mari

					}

					ftx = - kt * particles[n1].xix[n] - dampCoeffTan * v12x ;
					fty = - kt * particles[n1].xiy[n] - dampCoeffTan * v12y ;
					ftz = - kt * particles[n1].xiz[n] - dampCoeffTan * v12z ;



					fnn = sqrt(fnx*fnx + fny*fny + fnz*fnz);
					ftt = sqrt(ftx*ftx + fty*fty + ftz*ftz);

					fri = mu*fnn;

					t12x = ftx/ftt;
					t12y = fty/ftt;
					t12z = ftz/ftt;

					particles[n1].xix[n] += v12x*deltat;
					particles[n1].xiy[n] += v12y*deltat;
					particles[n1].xiz[n] += v12z*deltat;

					if (fri < ftt){

						ftx = fri*t12x;
						fty = fri*t12y;
						ftz = fri*t12z;

						particles[n1].xix[n] = - (1/kt) * (ftx + dampCoeffTan * v12x) ;
						particles[n1].xiy[n] = - (1/kt) * (fty + dampCoeffTan * v12y) ;
						particles[n1].xiz[n] = - (1/kt) * (ftz + dampCoeffTan * v12z) ;

					}

          felx += ftx;
          fely += fty;
          felz += ftz;

					}
					else {

						 fnx = sphereMat.linearStiff*del*n12x - sphereMat.dampCoeff * gn12 * n12x ;
						 fny = sphereMat.linearStiff*del*n12y - sphereMat.dampCoeff * gn12 * n12y ;
						 fnz = sphereMat.linearStiff*del*n12z - sphereMat.dampCoeff * gn12 * n12z ;

             felx = sphereMat.linearStiff*del*n12x;
             fely = sphereMat.linearStiff*del*n12y;
             felz = sphereMat.linearStiff*del*n12z;

             fvisx = - sphereMat.dampCoeff  * gn12 * n12x;
             fvisy = - sphereMat.dampCoeff  * gn12 * n12y;
             fvisz = - sphereMat.dampCoeff  * gn12 * n12z;

						 fnn = sqrt(fnx*fnx + fny*fny + fnz*fnz);


						 ftx = 0.0;
						 fty = 0.0;
						 ftz = 0.0;

					}


					int n_con = numContact;

					x1 = box(particles[n1].x0.x,demSize[0]);
					y1 = box(particles[n1].x0.y,demSize[1]);
					z1 = box(particles[n1].x0.z,demSize[2]);

					x2 = box(particles[n2].x0.x,demSize[0]);
					y2 = box(particles[n2].x0.y,demSize[1]);
					z2 = box(particles[n2].x0.z,demSize[2]);

					double zij = z1-z2;
					double xij = x1-x2;
					double yij = y1-y2;
					double rij = sqrt(zij*zij+yij*yij+xij*xij);

						if (rij <= 2.0*(maxPartRadius)){

							cPointIx[n_con] = x1;
							cPointIy[n_con] = y1;
							cPointIz[n_con] = z1;

							cPointJx[n_con] = x2;
							cPointJy[n_con] = y2;
							cPointJz[n_con] = z2;

							contactForceN[n_con] = fnn;

							birthTime[n_con] = (int) fnn;

							avgForce += fnn;

							numContact++;

              particles[n1].coord += 1.0;
              particles[n2].coord += 1.0;

						}


				}
				else {

					fnx = fny = fnz = ftx = fty = ftz = 0.0;
          felx = fely = felz = fvisx = fvisy = fvisz =0.0;

					particles[n1].xix[n] = 0.0;
					particles[n1].xiy[n] = 0.0;
					particles[n1].xiz[n] = 0.0;

				}

				dfx = fnx + ftx;
				dfy = fny + fty;
				dfz = fnz + ftz;

				dtx = fty*n12z - ftz*n12y;
				dty = ftz*n12x - ftx*n12z;
			  dtz = ftx*n12y - fty*n12x;

				particles[n1].FContact.x += dfx;
				particles[n1].FContact.y += dfy;
				particles[n1].FContact.z += dfz;

				particles[n2].FContact.x -= dfx;
				particles[n2].FContact.y -= dfy;
				particles[n2].FContact.z -= dfz;

        particles[n1].FElas = particles[n1].FElas + tVect(felx,fely,felz);
        particles[n1].FVisc = particles[n1].FVisc + tVect(fvisx,fvisy,fvisz);
        particles[n2].FElas = particles[n2].FElas - tVect(felx,fely,felz);
        particles[n2].FVisc = particles[n2].FVisc - tVect(fvisx,fvisy,fvisz);

				particles[n1].MContact.x += particles[n1].r * dtx ;
				particles[n1].MContact.y += particles[n1].r * dty ;
				particles[n1].MContact.z += particles[n1].r * dtz ;

				particles[n2].MContact.x += particles[n2].r * dtx ;
				particles[n2].MContact.y += particles[n2].r * dty ;
				particles[n2].MContact.z += particles[n2].r * dtz ;

				particles[n1].SContact.m00 += -x12  * dfx;
				particles[n1].SContact.m01 += -x12  * dfy;
				particles[n1].SContact.m02 += -x12  * dfz;
				particles[n1].SContact.m10 += -y12  * dfx;
				particles[n1].SContact.m11 += -y12  * dfy;
				particles[n1].SContact.m12 += -y12  * dfz;
				particles[n1].SContact.m20 += -z12  * dfx;
				particles[n1].SContact.m21 += -z12  * dfy;
				particles[n1].SContact.m22 += -z12  * dfz;

				particles[n2].SContact.m00 += -x12  * dfx;
				particles[n2].SContact.m01 += -x12  * dfy;
				particles[n2].SContact.m02 += -x12  * dfz;
				particles[n2].SContact.m10 += -y12  * dfx;
				particles[n2].SContact.m11 += -y12  * dfy;
				particles[n2].SContact.m12 += -y12  * dfz;
				particles[n2].SContact.m20 += -z12  * dfx;
				particles[n2].SContact.m21 += -z12  * dfy;
				particles[n2].SContact.m22 += -z12  * dfz;

        particles[n1].SElas.m00 += -x12  * felx;
        particles[n1].SElas.m01 += -x12  * fely;
        particles[n1].SElas.m02 += -x12  * felz;
        particles[n1].SElas.m10 += -y12  * felx;
        particles[n1].SElas.m11 += -y12  * fely;
        particles[n1].SElas.m12 += -y12  * felz;
        particles[n1].SElas.m20 += -z12  * felx;
        particles[n1].SElas.m21 += -z12  * fely;
        particles[n1].SElas.m22 += -z12  * felz;

        particles[n1].SVisc.m00 += -x12  * fvisx;
        particles[n1].SVisc.m01 += -x12  * fvisy;
        particles[n1].SVisc.m02 += -x12  * fvisz;
        particles[n1].SVisc.m10 += -y12  * fvisx;
        particles[n1].SVisc.m11 += -y12  * fvisy;
        particles[n1].SVisc.m12 += -y12  * fvisz;
        particles[n1].SVisc.m20 += -z12  * fvisx;
        particles[n1].SVisc.m21 += -z12  * fvisy;
        particles[n1].SVisc.m22 += -z12  * fvisz;

        particles[n2].SElas.m00 += -x12  * felx;
        particles[n2].SElas.m01 += -x12  * fely;
        particles[n2].SElas.m02 += -x12  * felz;
        particles[n2].SElas.m10 += -y12  * felx;
        particles[n2].SElas.m11 += -y12  * fely;
        particles[n2].SElas.m12 += -y12  * felz;
        particles[n2].SElas.m20 += -z12  * felx;
        particles[n2].SElas.m21 += -z12  * fely;
        particles[n2].SElas.m22 += -z12  * felz;

        particles[n2].SVisc.m00 += -x12  * fvisx;
        particles[n2].SVisc.m01 += -x12  * fvisy;
        particles[n2].SVisc.m02 += -x12  * fvisz;
        particles[n2].SVisc.m10 += -y12  * fvisx;
        particles[n2].SVisc.m11 += -y12  * fvisy;
        particles[n2].SVisc.m12 += -y12  * fvisz;
        particles[n2].SVisc.m20 += -z12  * fvisx;
        particles[n2].SVisc.m21 += -z12  * fvisy;
        particles[n2].SVisc.m22 += -z12  * fvisz;

				contactStress.m00 += -x12  * dfx;
				contactStress.m01 += -x12  * dfy;
				contactStress.m02 += -x12  * dfz;
				contactStress.m10 += -y12  * dfx;
				contactStress.m11 += -y12  * dfy;
				contactStress.m12 += -y12  * dfz;
				contactStress.m20 += -z12  * dfx;
				contactStress.m21 += -z12  * dfy;
				contactStress.m22 += -z12  * dfz;

			}
		}
	}

	avgForce /= (double) numContact ;

}

// Determining 2-simplicial complex (triangle) from contact network for persistent homology analysis by Perseus
void DEM::calculateSimplicialComplex(){

	double x1, x2, x3, y1, y2, y3, z1, z2, z3;

	numTriangle = 0;
//	numVertex = 0;

	for (int n1 = 0 ; n1 < particles.size(); n1++){

		for (int n = 1; n <= particles[n1].list[0] ; n++){

			int n2 = particles[n1].list[n];

			if (n2 < n1){

				double rad12 = 0.5*(particles[n1].r + particles[n2].r);
				double x12 = particles[n1].x0.x - particles[n2].x0.x ;
				double y12 = particles[n1].x0.y - particles[n2].x0.y ;
				double z12 = particles[n1].x0.z - particles[n2].x0.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

				double urx12 = particles[n1].x1.x - particles[n2].x1.x;
				double ury12 = particles[n1].x1.y - particles[n2].x1.y;
				double urz12 = particles[n1].x1.z - particles[n2].x1.z;

				double	n12x = x12/r12;
				double	n12y = y12/r12;
				double	n12z = z12/r12;


				double gn12 = urx12*n12x + ury12*n12y + urz12*n12z;

				double del12 = 2*rad12 - r12;

				if (del12 > 0.0) {

					for (int _n = 2; _n <= particles[n1].list[0] ; _n++){

					  int n3 = particles[n1].list[_n];

	 				  if (n3 < n2){

						  double rad13 = 0.5*(particles[n1].r + particles[n3].r);
						  double rad23 = 0.5*(particles[n2].r + particles[n3].r);

						  double x13 = particles[n1].x0.x - particles[n3].x0.x ;
						  double y13 = particles[n1].x0.y - particles[n3].x0.y ;
						  double z13 = particles[n1].x0.z - particles[n3].x0.z ;

						  double x23 = particles[n2].x0.x - particles[n3].x0.x ;
						  double y23 = particles[n2].x0.y - particles[n3].x0.y ;
						  double z23 = particles[n2].x0.z - particles[n3].x0.z ;

						  double r23 = sqrt(x23*x23 + y23*y23 + z23*z23);
						  double r13 = sqrt(x13*x13 + y13*y13 + z13*z13);

						  double urx13 = particles[n1].x1.x - particles[n3].x1.x;
						  double ury13 = particles[n1].x1.y - particles[n3].x1.y;
						  double urz13 = particles[n1].x1.z - particles[n3].x1.z;

						  double urx23 = particles[n2].x1.x - particles[n3].x1.x;
						  double ury23 = particles[n2].x1.y - particles[n3].x1.y;
						  double urz23 = particles[n2].x1.z - particles[n3].x1.z;

						  double	n13x = x13/r13;
						  double	n13y = y13/r13;
						  double	n13z = z13/r13;

						  double	n23x = x23/r23;
						  double	n23y = y23/r23;
						  double	n23z = z23/r23;

						  double gn23 = urx23*n23x + ury23*n23y + urz23*n23z;


						  double gn13 = urx13*n13x + ury13*n13y + urz13*n13z;

						  double del13 = 2*rad13 - r13;

						  double del23 = 2*rad23 - r23;


						  if (del13 > 0.0 && del23 > 0.0) {

							  int n_tri = numTriangle;

							  double fnx12 = sphereMat.linearStiff*del12*n12x - sphereMat.dampCoeff * gn12 * n12x ;
		 					  double fny12 = sphereMat.linearStiff*del12*n12y - sphereMat.dampCoeff * gn12 * n12y ;
		 					  double fnz12 = sphereMat.linearStiff*del12*n12z - sphereMat.dampCoeff * gn12 * n12z ;

							  double fnx13 = sphereMat.linearStiff*del13*n13x - sphereMat.dampCoeff * gn13 * n13x ;
		 					  double fny13 = sphereMat.linearStiff*del13*n13y - sphereMat.dampCoeff * gn13 * n13y ;
		 					  double fnz13 = sphereMat.linearStiff*del13*n13z - sphereMat.dampCoeff * gn13 * n13z ;

		 					  double fnx23 = sphereMat.linearStiff*del23*n23x - sphereMat.dampCoeff * gn23 * n23x ;
		 					  double fny23 = sphereMat.linearStiff*del23*n23y - sphereMat.dampCoeff * gn23 * n23y ;
		 					  double fnz23 = sphereMat.linearStiff*del23*n23z - sphereMat.dampCoeff * gn23 * n23z ;

							  double fnn12 = sqrt(fnx12*fnx12 + fny12*fny12 + fnz12*fnz12);

		 					  double fnn13 = sqrt(fnx13*fnx13 + fny13*fny13 + fnz13*fnz13);

		  					  double fnn23 = sqrt(fnx23*fnx23 + fny23*fny23 + fnz23*fnz23);

							  double forceMag = max(fnn12,max(fnn13,fnn23));

							  x1 = box(particles[n1].x0.x,demSize[0]);
							  y1 = box(particles[n1].x0.y,demSize[1]);
							  z1 = box(particles[n1].x0.z,demSize[2]);

							  x2 = box(particles[n2].x0.x,demSize[0]);
							  y2 = box(particles[n2].x0.y,demSize[1]);
							  z2 = box(particles[n2].x0.z,demSize[2]);

							  x3 = box(particles[n3].x0.x,demSize[0]);
							  y3 = box(particles[n3].x0.y,demSize[1]);
							  z3 = box(particles[n3].x0.z,demSize[2]);

							  vertexIx[n_tri] = x1;
							  vertexIy[n_tri] = y1;
							  vertexIz[n_tri] = z1;
						//	  numVertex++;
						//	  vertexI[n_tri] = numVertex;

							  vertexJx[n_tri] = x2;
							  vertexJy[n_tri] = y2;
							  vertexJz[n_tri] = z2;
						//	  numVertex++;
						//	  vertexJ[n_tri] = numVertex;

							  vertexKx[n_tri] = x3;
							  vertexKy[n_tri] = y3;
							  vertexKz[n_tri] = z3;
						//	  numVertex++;
						//	  vertexK[n_tri] = numVertex;

						//	  double fg = abs(particles[0].m * demF.y);

							  birthTime[n_tri] = (int) forceMag ;

							  numTriangle++;

						 	 }
	 					}
 					}
				}
			}
		}
	}
}

// Calculate contact forces on walls
void DEM::calculateWallContactForces() {
	double o1wx,o1wy,o1wz, v1wx, v1wy, v1wz, xi1w, fnxw, fnyw,fnzw, ktw, dampCoeffTanw, ftxw, ftyw, ftzw;
	double fnnw, fttw, friw, t1wx, t1wy, t1wz, dfxw, dfyw, dfzw, dtxw, dtyw, dtzw;
  double felx, fely, felz, fvisx, fvisy, fvisz;
	double n1wx, n1wy, n1wz, gn1w, muw;

	for (int n1 = 0; n1 < particles.size(); n1++){

		for (int n2_w = 0; n2_w < walls.size(); n2_w++){

			double radw = particles[n1].r ;
			double x1w = (particles[n1].x0.x - walls[n2_w].p.x)*walls[n2_w].n.x*walls[n2_w].n.x ;
			double y1w = (particles[n1].x0.y - walls[n2_w].p.y)*walls[n2_w].n.y*walls[n2_w].n.y ;
			double z1w = (particles[n1].x0.z - walls[n2_w].p.z)*walls[n2_w].n.z*walls[n2_w].n.z ;
			double r1w = sqrt(x1w*x1w + y1w*y1w + z1w*z1w);

			double urxw = particles[n1].x1.x ;
			double uryw = particles[n1].x1.y ;
			double urzw = particles[n1].x1.z ;

			double delw = radw - r1w;
		//	double drw = -delw;


			if (delw > 0.0){

				n1wx = x1w/r1w;
			    n1wy = y1w/r1w;
				n1wz = z1w/r1w;


				 gn1w = urxw*n1wx + uryw*n1wy + urzw*n1wz;

				 muw = sphereMat.frictionCoefPart;

			if (muw != 0.0 ){

				o1wx = particles[n1].r * particles[n1].w1.x;
				o1wy = particles[n1].r * particles[n1].w1.y;
				o1wz = particles[n1].r * particles[n1].w1.z;

				v1wx = urxw + n1wy*o1wz - n1wz*o1wy - gn1w*n1wx;
				v1wy = uryw + n1wz*o1wx - n1wx*o1wz - gn1w*n1wy;
				v1wz = urzw + n1wx*o1wy - n1wy*o1wx - gn1w*n1wz;

				xi1w = particles[n1].xiw.x*n1wx + particles[n1].xiw.y*n1wy + particles[n1].xiw.z*n1wz ;

				particles[n1].xiw.x -= xi1w*n1wx;
				particles[n1].xiw.y -= xi1w*n1wy;
				particles[n1].xiw.z -= xi1w*n1wz;

				fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
				fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
				fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;

        felx = sphereMat.linearStiff*delw*n1wx;
        fely = sphereMat.linearStiff*delw*n1wy;
        felz = sphereMat.linearStiff*delw*n1wz;

        fvisx = - sphereMat.dampCoeff * gn1w * n1wx;
        fvisy = - sphereMat.dampCoeff * gn1w * n1wy;
        fvisz = - sphereMat.dampCoeff * gn1w * n1wz;


				//tangential stiffness
				ktw = 0.2*sphereMat.linearStiff;


				// no tangential dissipation as Seto et al
				if (sphereMat.kinVisc == 0.0){
				dampCoeffTanw = sphereMat.dampCoeff;
				}
				else {
				dampCoeffTanw = 0.0;
				}

				ftxw = - ktw * particles[n1].xiw.x - dampCoeffTanw * v1wx ;
				ftyw = - ktw * particles[n1].xiw.y - dampCoeffTanw * v1wy ;
				ftzw = - ktw * particles[n1].xiw.z - dampCoeffTanw * v1wz ;

				fnnw = sqrt(fnxw*fnxw + fnyw*fnyw + fnzw*fnzw);
				fttw = sqrt(ftxw*ftxw + ftyw*ftyw + ftzw*ftzw);

				friw = muw*fnnw;

				t1wx = ftxw/fttw;
				t1wy = ftyw/fttw;
				t1wz = ftzw/fttw;

				particles[n1].xiw.x += v1wx*deltat;
				particles[n1].xiw.y += v1wy*deltat;
				particles[n1].xiw.z += v1wz*deltat;

				if (friw < fttw){

					ftxw = friw*t1wx;
					ftyw = friw*t1wy;
					ftzw = friw*t1wz;

					particles[n1].xiw.x = - (1/ktw) * (ftxw + dampCoeffTanw * v1wx) ;
					particles[n1].xiw.y = - (1/ktw) * (ftyw + dampCoeffTanw * v1wy) ;
					particles[n1].xiw.z = - (1/ktw) * (ftzw + dampCoeffTanw * v1wz) ;

				}
        felx += ftxw;
        fely += ftyw;
        felz += ftzw;

				}
				else {

					 fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
					 fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
					 fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;

           felx = sphereMat.linearStiff*delw*n1wx;
           fely = sphereMat.linearStiff*delw*n1wy;
           felz = sphereMat.linearStiff*delw*n1wz;

           fvisx = - sphereMat.dampCoeff * gn1w * n1wx;
           fvisy = - sphereMat.dampCoeff * gn1w * n1wy;
           fvisz = - sphereMat.dampCoeff * gn1w * n1wz;

					 ftxw = 0.0;
					 ftyw = 0.0;
					 ftzw = 0.0;

				}
			}
			else {

				fnxw = fnyw = fnzw = ftxw = ftyw = ftzw = 0.0;
        felx = fely = felz = fvisx = fvisy = fvisz =0.0;

				particles[n1].xiw.x = 0.0;
				particles[n1].xiw.y = 0.0;
				particles[n1].xiw.z = 0.0;

			}

			double dfxw = fnxw + ftxw;
			double dfyw = fnyw + ftyw;
			double dfzw = fnzw + ftzw;


			double dtxw = ftyw*n1wz - ftzw*n1wy;
			double dtyw = ftzw*n1wx - ftxw*n1wz;
			double dtzw = ftxw*n1wy - ftyw*n1wx;

			particles[n1].FContact.x += dfxw;
			particles[n1].FContact.y += dfyw;
			particles[n1].FContact.z += dfzw;

			walls[n2_w].FContact.x -= dfxw;
			walls[n2_w].FContact.y -= dfyw;
			walls[n2_w].FContact.z -= dfzw;

      particles[n1].FElas = particles[n1].FElas + tVect(felx,fely,felz);
      particles[n1].FVisc = particles[n1].FVisc + tVect(fvisx,fvisy,fvisz);

			particles[n1].MContact.x += particles[n1].r * dtxw ;
			particles[n1].MContact.y += particles[n1].r * dtyw ;
			particles[n1].MContact.z += particles[n1].r * dtzw ;

			particles[n1].SContact.m00 += -x1w  * dfxw;
			particles[n1].SContact.m01 += -x1w  * dfyw;
			particles[n1].SContact.m02 += -x1w  * dfzw;
			particles[n1].SContact.m10 += -y1w  * dfxw;
			particles[n1].SContact.m11 += -y1w  * dfyw;
			particles[n1].SContact.m12 += -y1w  * dfzw;
			particles[n1].SContact.m20 += -z1w  * dfxw;
			particles[n1].SContact.m21 += -z1w  * dfyw;
			particles[n1].SContact.m22 += -z1w  * dfzw;

      particles[n1].SElas.m00 += -x1w  * felx;
      particles[n1].SElas.m01 += -x1w  * fely;
      particles[n1].SElas.m02 += -x1w  * felz;
      particles[n1].SElas.m10 += -y1w  * felx;
      particles[n1].SElas.m11 += -y1w  * fely;
      particles[n1].SElas.m12 += -y1w  * felz;
      particles[n1].SElas.m20 += -z1w  * felx;
      particles[n1].SElas.m21 += -z1w  * fely;
      particles[n1].SElas.m22 += -z1w  * felz;

      particles[n1].SVisc.m00 += -x1w  * fvisx;
      particles[n1].SVisc.m01 += -x1w  * fvisy;
      particles[n1].SVisc.m02 += -x1w  * fvisz;
      particles[n1].SVisc.m10 += -y1w  * fvisx;
      particles[n1].SVisc.m11 += -y1w  * fvisy;
      particles[n1].SVisc.m12 += -y1w  * fvisz;
      particles[n1].SVisc.m20 += -z1w  * fvisx;
      particles[n1].SVisc.m21 += -z1w  * fvisy;
      particles[n1].SVisc.m22 += -z1w  * fvisz;

			contactStress.m00 += -x1w  * dfxw;
			contactStress.m01 += -x1w  * dfyw;
			contactStress.m02 += -x1w  * dfzw;
			contactStress.m10 += -y1w  * dfxw;
			contactStress.m11 += -y1w  * dfyw;
			contactStress.m12 += -y1w  * dfzw;
			contactStress.m20 += -z1w  * dfxw;
			contactStress.m21 += -z1w  * dfyw;
			contactStress.m22 += -z1w  * dfzw;

		}
	}

}



void DEM::calculateFictionalWallContactForces() {
	double o1wx,o1wy,o1wz, v1wx, v1wy, v1wz, xi1w, fnxw, fnyw,fnzw, ktw, dampCoeffTanw, ftxw, ftyw, ftzw;
	double fnnw, fttw, friw, t1wx, t1wy, t1wz, dfxw, dfyw, dfzw, dtxw, dtyw, dtzw;

	double n1wx, n1wy, n1wz, gn1w, muw;

	for (int n1 = 0; n1 < particles.size(); n1++){


			double radw = particles[n1].r ;
			double x1w = 0.0 ;
      double y1w = particles[n1].x0.y - (demSize[1] - 0.3);
			double z1w = 0.0;
			double r1w = sqrt(x1w*x1w + y1w*y1w + z1w*z1w);

			double urxw = particles[n1].x1.x ;
			double uryw = particles[n1].x1.y ;
			double urzw = particles[n1].x1.z ;

			double delw = radw - r1w;
		//	double drw = -delw;


			if (delw > 0.0){

				n1wx = x1w/r1w;
			  n1wy = y1w/r1w;
				n1wz = z1w/r1w;


				 gn1w = urxw*n1wx + uryw*n1wy + urzw*n1wz;

				 muw = sphereMat.frictionCoefPart;

			if (muw != 0.0 ){

				o1wx = particles[n1].r * particles[n1].w1.x;
				o1wy = particles[n1].r * particles[n1].w1.y;
				o1wz = particles[n1].r * particles[n1].w1.z;

				v1wx = urxw + n1wy*o1wz - n1wz*o1wy - gn1w*n1wx;
				v1wy = uryw + n1wz*o1wx - n1wx*o1wz - gn1w*n1wy;
				v1wz = urzw + n1wx*o1wy - n1wy*o1wx - gn1w*n1wz;

				xi1w = particles[n1].xiw.x*n1wx + particles[n1].xiw.y*n1wy + particles[n1].xiw.z*n1wz ;

				particles[n1].xiw.x -= xi1w*n1wx;
				particles[n1].xiw.y -= xi1w*n1wy;
				particles[n1].xiw.z -= xi1w*n1wz;

				fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
				fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
				fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;

				//tangential stiffness
				ktw = 0.2*sphereMat.linearStiff;


				// no tangential dissipation as Seto et al
				if (sphereMat.kinVisc == 0.0){
				dampCoeffTanw = sphereMat.dampCoeff;
				}
				else {
				dampCoeffTanw = 0.0;
				}

				ftxw = - ktw * particles[n1].xiw.x - dampCoeffTanw * v1wx ;
				ftyw = - ktw * particles[n1].xiw.y - dampCoeffTanw * v1wy ;
				ftzw = - ktw * particles[n1].xiw.z - dampCoeffTanw * v1wz ;

				fnnw = sqrt(fnxw*fnxw + fnyw*fnyw + fnzw*fnzw);
				fttw = sqrt(ftxw*ftxw + ftyw*ftyw + ftzw*ftzw);

				friw = muw*fnnw;

				t1wx = ftxw/fttw;
				t1wy = ftyw/fttw;
				t1wz = ftzw/fttw;

				particles[n1].xiw.x += v1wx*deltat;
				particles[n1].xiw.y += v1wy*deltat;
				particles[n1].xiw.z += v1wz*deltat;

				if (friw < fttw){

					ftxw = friw*t1wx;
					ftyw = friw*t1wy;
					ftzw = friw*t1wz;

					particles[n1].xiw.x = - (1/ktw) * (ftxw + dampCoeffTanw * v1wx) ;
					particles[n1].xiw.y = - (1/ktw) * (ftyw + dampCoeffTanw * v1wy) ;
					particles[n1].xiw.z = - (1/ktw) * (ftzw + dampCoeffTanw * v1wz) ;

				}


				}
				else {

					 fnxw = sphereMat.linearStiff*delw*n1wx - sphereMat.dampCoeff * gn1w * n1wx ;
					 fnyw = sphereMat.linearStiff*delw*n1wy - sphereMat.dampCoeff * gn1w * n1wy ;
					 fnzw = sphereMat.linearStiff*delw*n1wz - sphereMat.dampCoeff * gn1w * n1wz ;


					 ftxw = 0.0;
					 ftyw = 0.0;
					 ftzw = 0.0;

				}
			}
			else {

				fnxw = fnyw = fnzw = ftxw = ftyw = ftzw = 0.0;

				particles[n1].xiw.x = 0.0;
				particles[n1].xiw.y = 0.0;
				particles[n1].xiw.z = 0.0;

			}

			double dfxw = fnxw + ftxw;
			double dfyw = fnyw + ftyw;
			double dfzw = fnzw + ftzw;


			double dtxw = ftyw*n1wz - ftzw*n1wy;
			double dtyw = ftzw*n1wx - ftxw*n1wz;
			double dtzw = ftxw*n1wy - ftyw*n1wx;

			particles[n1].FContact.x += dfxw;
			particles[n1].FContact.y += dfyw;
			particles[n1].FContact.z += dfzw;

			particles[n1].MContact.x += particles[n1].r * dtxw ;
			particles[n1].MContact.y += particles[n1].r * dtyw ;
			particles[n1].MContact.z += particles[n1].r * dtzw ;

	}

}




// Calculate lubrication forces on particles from Nguyen and Ladd, Kim and Karilla 1991
void DEM::calculateParticleLubricationForces() {
  // Also from Susp3D, based on Nguyen and Ladd, PRE (2002)
	double        Xa, Ya, Yb1, Yb2, Xg1, Xg2, Yg1, Yg2;
	double        Yc11, Yc12, Yc21, Yc22, Yh11, Yh12, Yh21, Yh22;
	double        Fx, Fy, Fz, T1x, T1y, T1z, T2x, T2y, T2z;
	double        S1xx, S1yy, S1zz, S1yz, S1zx, S1xy;
	double        S2xx, S2yy, S2zz, S2yz, S2zx, S2xy;
	double        nu, rad, cut, dr, h, h_c;
	double        umx, umy, umz, w1x, w1y, w1z, w2x, w2y, w2z;
	double        umxrx, umxry, umxrz, w1xrx, w1xry, w1xrz, w2xrx, w2xry, w2xrz;
	double        um_r, w1_r, w2_r, beta;

	double lub_N = sphereMat.lubNormal;
	double lub_T = sphereMat.lubTangential;
	double lub_R = sphereMat.lubRotation;


	for (int n1 =0; n1 < particles.size(); n1++){

		for (int n = 1; n <= particles[n1].list[0]; n++){

			int n2 = particles[n1].list[n];

			if (n2 < n1){

				cut = sphereMat.lubCutoff;
				rad = 0.5*(particles[n1].r + particles[n2].r);
				double x12 = particles[n1].x0.x - particles[n2].x0.x ;
				double y12 = particles[n1].x0.y - particles[n2].x0.y ;
				double z12 = particles[n1].x0.z - particles[n2].x0.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

				umx = particles[n1].x1.x - particles[n2].x1.x;
				umy = particles[n1].x1.y - particles[n2].x1.y;
				umz = particles[n1].x1.z - particles[n2].x1.z;

				w1x   = particles[n1].w1.x;
				w1y   = particles[n1].w1.y;
				w1z   = particles[n1].w1.z;
				w2x   = particles[n2].w1.x;
				w2y   = particles[n2].w1.y;
				w2z   = particles[n2].w1.z;


				double dr = r12 - 2.0*rad;
				double h = max(dr,cut);

				if ( dr < max(cut, max(lub_N, max(lub_T,lub_R))))
				{


					double nu = sphereMat.kinVisc;
					double alpha = M_PI*nu*rad*deltat;
					double beta1 = particles[n2].r / particles[n1].r;
					double beta2 = beta1*beta1;
					double beta3 = beta1*beta2;
					double beta4 = pow((1+beta1), 4.0);
					double beta5 = pow((1+beta1), 5.0);
					if (dr < max(cut, lub_N))
					{
            particles[n1].lub_coord += 1.0;
            particles[n2].lub_coord += 1.0;
						Xa   =  4.0*beta2/beta4;                         /* Set lubrication coefficients */
						Xg1  = 12.0*beta2/beta5;
						Xg2  =-12.0*beta3/beta5;
						h_c   = max(dr, lub_N);
						Xa   *= rad/h - rad/h_c;
						Xg1  *= rad/h - rad/h_c;
						Xg2  *= rad/h - rad/h_c;
						Xa   *= 6.0*alpha;
						Xg1  *= 4.0*alpha*rad;
						Xg2  *= 4.0*alpha*rad;
					}
					else
						Xa = Xg1 = Xg2 = 0.0;
					if (dr < max(cut, lub_T))
					{
						Ya    = 8.0*beta1*(2.0+beta1+2.0*beta2)/(15.0*beta4);
						Yb1   =-4.0*beta1*(4.0+beta1)/(5.0*beta4);
						Yb2   = 4.0*beta1*(4.0*beta2+beta1)/(5.0*beta4);
						Yg1   = 2.0*beta1*(4.0-beta1+7.0*beta2)/(5.0*beta5);
						Yg2   =-2.0*beta2*(4.0*beta2-beta1+7.0)/(5.0*beta5);
						h_c   = max(dr, lub_T);
						beta  = alpha*log(h_c/h);
						Ya   *= 6.0*beta;
						Yb1  *= 4.0*beta*rad;
						Yb2  *= 4.0*beta*rad;
						Yg1  *= 4.0*beta*rad;
						Yg2  *= 4.0*beta*rad;
					}
					else
						Ya = Yb1 = Yb2 = Yg1 = Yg2 = 0.0;
					if (dr < max(cut, lub_R))
					{
						Yc11  = 16.0*beta1/(5.0*beta4);
						Yc12  =  4.0*beta2/(5.0*beta4);
						Yc21  =  4.0*beta2/(5.0*beta4);
						Yc22  = 16.0*beta3/(5.0*beta4);
						Yh11  =  4.0*beta1*(2.0-beta1)/(5.0*beta5);
						Yh12  =  2.0*beta2*(1.0+7.0*beta1)/(5.0*beta5);
						Yh21  =  2.0*beta2*(beta1+7.0)/(5.0*beta5);
						Yh22  =  4.0*beta2*(2.0*beta2-beta1)/(5.0*beta5);
						h_c   = max(dr, lub_R);
						beta  = alpha*log(h_c/h);
						Yc11 *= 8.0*beta*rad*rad;
						Yc12 *= 8.0*beta*rad*rad;
						Yc21 *= 8.0*beta*rad*rad;
						Yc22 *= 8.0*beta*rad*rad;
						Yh11 *= 8.0*beta*rad*rad;
						Yh21 *= 8.0*beta*rad*rad;
						Yh12 *= 8.0*beta*rad*rad;
						Yh22 *= 8.0*beta*rad*rad;
					}
					else
						Yc11 = Yc12 = Yc21 = Yc22 = Yh11 = Yh12 = Yh21 = Yh22 = 0.0;

					x12  /= r12;
					y12  /= r12;
					z12  /= r12;


					um_r  = umx*x12 + umy*y12 + umz*z12;
					w1_r  = w1x*x12 + w1y*y12 + w1z*z12;
					w2_r  = w2x*x12 + w2y*y12 + w2z*z12;
					umxrx = umy*z12 - umz*y12;
					umxry = umz*x12 - umx*z12;
					umxrz = umx*y12 - umy*x12;
					w1xrx = w1y*z12 - w1z*y12;
					w1xry = w1z*x12 - w1x*z12;
					w1xrz = w1x*y12 - w1y*x12;
					w2xrx = w2y*z12 - w2z*y12;
					w2xry = w2z*x12 - w2x*z12;
					w2xrz = w2x*y12 - w2y*x12;

					Fx    = -Xa*um_r*x12 - Ya*(umx-um_r*x12);
					Fy    = -Xa*um_r*y12 - Ya*(umy-um_r*y12);
					Fz    = -Xa*um_r*z12 - Ya*(umz-um_r*z12);
					Fx   += -Yb1*w1xrx + Yb2*w2xrx;
					Fy   += -Yb1*w1xry + Yb2*w2xry;
					Fz   += -Yb1*w1xrz + Yb2*w2xrz;
					T1x   =  Yb1*umxrx;
					T1y   =  Yb1*umxry;
					T1z   =  Yb1*umxrz;
					T2x   = -Yb2*umxrx;
					T2y   = -Yb2*umxry;
					T2z   = -Yb2*umxrz;
					T1x  += -Yc11*(w1x-w1_r*x12) - Yc12*(w2x-w2_r*x12);
					T1y  += -Yc11*(w1y-w1_r*y12) - Yc12*(w2y-w2_r*y12);
					T1z  += -Yc11*(w1z-w1_r*z12) - Yc12*(w2z-w2_r*z12);
					T2x  += -Yc21*(w1x-w1_r*x12) - Yc22*(w2x-w2_r*x12);
					T2y  += -Yc21*(w1y-w1_r*y12) - Yc22*(w2y-w2_r*y12);
					T2z  += -Yc21*(w1z-w1_r*z12) - Yc22*(w2z-w2_r*z12);
					S1xx  =  Xg1*um_r*x12*x12 + Yg1*(x12*umx+umx*x12);
					S1yy  =  Xg1*um_r*y12*y12 + Yg1*(y12*umy+umy*y12);
					S1zz  =  Xg1*um_r*z12*z12 + Yg1*(z12*umz+umz*z12);
					S1yz  =  Xg1*um_r*y12*z12 + Yg1*(y12*umz+umy*z12);
					S1zx  =  Xg1*um_r*z12*x12 + Yg1*(z12*umx+umz*x12);
					S1xy  =  Xg1*um_r*x12*y12 + Yg1*(x12*umy+umx*y12);
					S2xx  = -Xg2*um_r*x12*x12 - Yg2*(x12*umx+umx*x12);
					S2yy  = -Xg2*um_r*y12*y12 - Yg2*(y12*umy+umy*y12);
					S2zz  = -Xg2*um_r*z12*z12 - Yg2*(z12*umz+umz*z12);
					S2yz  = -Xg2*um_r*y12*z12 - Yg2*(y12*umz+umy*z12);
					S2zx  = -Xg2*um_r*z12*x12 - Yg2*(z12*umx+umz*x12);
					S2xy  = -Xg2*um_r*x12*y12 - Yg2*(x12*umy+umx*y12);
					S1xx += -Yh11*(x12*w1xrx+w1xrx*x12) - Yh12*(x12*w2xrx+w2xrx*x12);
					S1yy += -Yh11*(y12*w1xry+w1xry*y12) - Yh12*(y12*w2xry+w2xry*y12);
					S1zz += -Yh11*(z12*w1xrz+w1xrz*z12) - Yh12*(z12*w2xrz+w2xrz*z12);
					S1yz += -Yh11*(y12*w1xrz+w1xry*z12) - Yh12*(y12*w2xrz+w2xry*z12);
					S1zx += -Yh11*(z12*w1xrx+w1xrz*x12) - Yh12*(z12*w2xrx+w2xrz*x12);
					S1xy += -Yh11*(x12*w1xry+w1xrx*y12) - Yh12*(x12*w2xry+w2xrx*y12);
					S2xx += -Yh21*(x12*w1xrx+w1xrx*x12) - Yh22*(x12*w2xrx+w2xrx*x12);
					S2yy += -Yh21*(y12*w1xry+w1xry*y12) - Yh22*(y12*w2xry+w2xry*y12);
					S2zz += -Yh21*(z12*w1xrz+w1xrz*z12) - Yh22*(z12*w2xrz+w2xrz*z12);
					S2yz += -Yh21*(y12*w1xrz+w1xry*z12) - Yh22*(y12*w2xrz+w2xry*z12);
					S2zx += -Yh21*(z12*w1xrx+w1xrz*x12) - Yh22*(z12*w2xrx+w2xrz*x12);
					S2xy += -Yh21*(x12*w1xry+w1xrx*y12) - Yh22*(x12*w2xry+w2xrx*y12);


					particles[n1].FLub.x += Fx;
					particles[n1].FLub.y += Fy;
					particles[n1].FLub.z += Fz;
					particles[n2].FLub.x -= Fx;
					particles[n2].FLub.y -= Fy;
					particles[n2].FLub.z -= Fz;
					particles[n1].MLub.x += T1x;
					particles[n1].MLub.y += T1y;
					particles[n1].MLub.z += T1z;
					particles[n2].MLub.x += T2x;
					particles[n2].MLub.y += T2y;
					particles[n2].MLub.z += T2z;


					particles[n1].SLub.m00 += S1xx;
					particles[n1].SLub.m01 += S1xy;
					particles[n1].SLub.m02 += S1zx;
					particles[n1].SLub.m10 += S1xy;
					particles[n1].SLub.m11 += S1yy;
					particles[n1].SLub.m12 += S1yz;
					particles[n1].SLub.m20 += S1zx;
					particles[n1].SLub.m21 += S1yz;
					particles[n1].SLub.m22 += S1zz;

					particles[n2].SLub.m00 += S2xx;
					particles[n2].SLub.m01 += S2xy;
					particles[n2].SLub.m02 += S2zx;
					particles[n2].SLub.m10 += S2xy;
					particles[n2].SLub.m11 += S2yy;
					particles[n2].SLub.m12 += S2yz;
					particles[n2].SLub.m20 += S2zx;
					particles[n2].SLub.m21 += S2yz;
					particles[n2].SLub.m22 += S2zz;



				}
			}
		}
	}

}


// Calculate lubrication forces on particles due to walls
void DEM::calculateWallLubricationForces()  {
	double        Xa, Ya, Yb1, Yb2, Xg1, Xg2, Yg1, Yg2;
 	double        Yc11, Yc12, Yc21, Yc22, Yh11, Yh12, Yh21, Yh22;
 	double        Fx, Fy, Fz, T1x, T1y, T1z, T2x, T2y, T2z;
 	double        S1xx, S1yy, S1zz, S1yz, S1zx, S1xy;
 	double        S2xx, S2yy, S2zz, S2yz, S2zx, S2xy;
 	double        nu, rad, cut, dr, h, h_c;
 	double        umx, umy, umz, w1x, w1y, w1z, w2x, w2y, w2z;
 	double        umxrx, umxry, umxrz, w1xrx, w1xry, w1xrz, w2xrx, w2xry, w2xrz;
 	double        um_r, w1_r, w2_r, beta;

 	double lub_N = sphereMat.lubNormal;
 	double lub_T = sphereMat.lubTangential;
 	double lub_R = sphereMat.lubRotation;


 	for (int n1 = 0; n1 < particles.size(); n1++){

			for (int n2_w = 0; n2_w < walls.size(); n2_w++){


			  cut = sphereMat.lubCutoff;
			  rad = particles[n1].r ;
				double x12 = (particles[n1].x0.x - walls[n2_w].p.x)*walls[n2_w].n.x*walls[n2_w].n.x ;
				double y12 = (particles[n1].x0.y - walls[n2_w].p.y)*walls[n2_w].n.y*walls[n2_w].n.y ;
				double z12 = (particles[n1].x0.z - walls[n2_w].p.z)*walls[n2_w].n.z*walls[n2_w].n.z ;
				double r12 = sqrt(x12*x12 + y12*y12 + z12*z12);

			 	umx = particles[n1].x1.x ;
			 	umy = particles[n1].x1.y ;
				umz = particles[n1].x1.z ;

				w1x   = particles[n1].w1.x;
				w1y   = particles[n1].w1.y;
				w1z   = particles[n1].w1.z;

				//r12 = sqrt(x12*x12 + y12*y12 + z12*z12);
				dr  = r12-rad;
				h   = max(dr, cut);


				if (dr < max(cut, max(lub_N, max(lub_T, lub_R))))
				{

					double nu = sphereMat.kinVisc;
					double alpha = M_PI*nu*rad*deltat;

					if (dr < max(cut, lub_N))
					{
						Xa    = 1.0;                                             /* Set lubrication coefficients */
						Xg1   = 3.0/2.0;
						h_c   = max(dr, lub_N);
						Xa   *= rad/h - rad/h_c;
						Xg1  *= rad/h - rad/h_c;
						Xa   *= 6.0*alpha;
						Xg1  *= 4.0*alpha*rad;
					}
					else Xa = Xg1 = 0.0;
					if (dr < max(cut, lub_T))
					{
						Ya    = 8.0/15.0;
						Yb1   =-1.0/5.0;
						Yg1   = 7.0/10.0;
						h_c   = max(dr, lub_T);
						beta  = alpha*log(h_c/h);
						Ya   *= 6.0*beta;
						Yb1  *= 4.0*beta*rad;
						Yg1  *= 4.0*beta*rad;
					}
					else Ya = Yb1 = Yg1 = 0.0;

					if (dr < max(cut, lub_R))
					{
						Yc11  = 2.0/5.0;
						Yh11  =-1.0/10.0;
						h_c   = max(dr, lub_R);
						beta  = alpha*log(h_c/h);
						Yc11 *= 8.0*beta*rad*rad;
						Yh11 *= 8.0*beta*rad*rad;
					}
					else Yc11 = Yh11 = 0.0;

					x12 /= r12;
					y12 /= r12;
					z12 /= r12;

					um_r = umx*x12 + umy*y12 + umz*z12;
					w1_r = w1x*x12 + w1y*y12 + w1z*z12;
					umxrx = umy*z12 - umz*y12;
					umxry = umz*x12 - umx*z12;
					umxrz = umx*y12 - umy*x12;
					w1xrx = w1y*z12 - w1z*y12;
					w1xry = w1z*x12 - w1x*z12;
					w1xrz = w1x*y12 - w1y*x12;

					Fx    = -Xa*x12*um_r - Ya*(umx-um_r*x12);
					Fy    = -Xa*y12*um_r - Ya*(umy-um_r*y12);
					Fz    = -Xa*z12*um_r - Ya*(umz-um_r*z12);
					Fx   += -Yb1*w1xrx;
					Fy   += -Yb1*w1xry;
					Fz   += -Yb1*w1xrz;
					T1x   =  Yb1*umxrx;
					T1y   =  Yb1*umxry;
					T1z   =  Yb1*umxrz;
					T1x  += -Yc11*(w1x-w1_r*x12);
					T1y  += -Yc11*(w1y-w1_r*y12);
					T1z  += -Yc11*(w1z-w1_r*z12);

					S1xx  =  Xg1*(um_r*x12*x12) + Yg1*(x12*umx+umx*x12);
					S1yy  =  Xg1*(um_r*y12*y12) + Yg1*(y12*umy+umy*y12);
					S1zz  =  Xg1*(um_r*z12*z12) + Yg1*(z12*umz+umz*z12);
					S1yz  =  Xg1*(um_r*y12*z12) + Yg1*(y12*umz+umy*z12);
					S1zx  =  Xg1*(um_r*z12*x12) + Yg1*(z12*umx+umz*x12);
					S1xy  =  Xg1*(um_r*x12*y12) + Yg1*(x12*umy+umx*y12);
					S1xx += -Yh11*(x12*w1xrx+w1xrx*x12);
					S1yy += -Yh11*(y12*w1xry+w1xry*y12);
					S1zz += -Yh11*(z12*w1xrz+w1xrz*z12);
					S1yz += -Yh11*(y12*w1xrz+w1xry*z12);
					S1zx += -Yh11*(z12*w1xrx+w1xrz*x12);
					S1xy += -Yh11*(x12*w1xry+w1xrx*y12);

					particles[n1].FLub.x += Fx;
					particles[n1].FLub.y += Fy;
					particles[n1].FLub.z += Fz;

					particles[n1].MLub.x += T1x;
					particles[n1].MLub.y += T1y;
					particles[n1].MLub.z += T1z;

					walls[n2_w].FLub.x -= Fx;
					walls[n2_w].FLub.y -= Fy;
					walls[n2_w].FLub.z -= Fz;

					particles[n1].SLub.m00 += S1xx;
					particles[n1].SLub.m01 += S1xy;
					particles[n1].SLub.m02 += S1zx;
					particles[n1].SLub.m10 += S1xy;
					particles[n1].SLub.m11 += S1yy;
					particles[n1].SLub.m12 += S1yz;
					particles[n1].SLub.m20 += S1zx;
					particles[n1].SLub.m21 += S1yz;
					particles[n1].SLub.m22 += S1zz;

			}

		}
 	}
}

// integration functions

void DEM::integration(int timeStep){

	double xi,yi,zi;
	double xo,yo,zo;

	for (int n =0 ; n < particles.size() ; n++){

		xi = particles[n].x0.x;
		yi = particles[n].x0.y;
		zi = particles[n].x0.z;


    // update position
    particles[n].x0.x += deltat * particles[n].x1.x + deltat * deltat * particles[n].FTotal.x * particles[n]._m;
    particles[n].x0.y += deltat * particles[n].x1.y + deltat * deltat * particles[n].FTotal.y * particles[n]._m;
    particles[n].x0.z += deltat * particles[n].x1.z + deltat * deltat * particles[n].FTotal.z * particles[n]._m;

    particles[n].x1.x += deltat * particles[n].FTotal.x * particles[n]._m;
    particles[n].x1.y += deltat * particles[n].FTotal.y * particles[n]._m;
    particles[n].x1.z += deltat * particles[n].FTotal.z * particles[n]._m;
		// update angular velocity
		particles[n].w1.x += deltat * particles[n].MTotal.x * particles[n]._I;
		particles[n].w1.y += deltat * particles[n].MTotal.y * particles[n]._I;
		particles[n].w1.z += deltat * particles[n].MTotal.z * particles[n]._I;

		particles[n]._x0.x = xi;
		particles[n]._x0.y = yi;
		particles[n]._x0.z = zi;

		particles[n].disp =  particles[n].x_i - particles[n].x0;

	}


	for (int n =0 ; n < particles.size() ; n++){

		particles[n].disp_r = sqrt(particles[n].disp.x*particles[n].disp.x + particles[n].disp.z*particles[n].disp.z);

	}

	for (int n =0 ; n < particles.size() ; n++){
		particles[n].trueStrain = log(particles[n].x0.y / particles[n]._x0.y );
		particles[n].engStrain = particles[n].disp.y / particles[n]._x0.y ;
	}

  if (timeStep == objWaitTime){
    for (int n =0 ; n < particles.size() ; n++){
    particles[n].x_i = particles[n].x0;
  }
  }

	if (timeStep >= objWaitTime) {


    for (int b = 0; b < bodies.size() ; b++){

      bodies[b].x0.x += deltat * bodies[b].x1.x + deltat * deltat * bodies[b].FTotal.x * bodies[b]._m;
      bodies[b].x0.y += deltat * bodies[b].x1.y + deltat * deltat * bodies[b].FTotal.y * bodies[b]._m;
      bodies[b].x0.z += deltat * bodies[b].x1.z + deltat * deltat * bodies[b].FTotal.z * bodies[b]._m;

      bodies[b].x1.x += deltat * bodies[b].FTotal.x * bodies[b]._m;
      bodies[b].x1.y += deltat * bodies[b].FTotal.y * bodies[b]._m;
      bodies[b].x1.z += deltat * bodies[b].FTotal.z * bodies[b]._m;

    }

    for (int o = 0 ; o < objects.size() ; o++){



      objects[o].x0.x += deltat * objects[o].x1.x + deltat * deltat * objects[o].FTotal.x * objects[o]._m;
			objects[o].x0.y += deltat * objects[o].x1.y + deltat * deltat * objects[o].FTotal.y * objects[o]._m;
			objects[o].x0.z += deltat * objects[o].x1.z + deltat * deltat * objects[o].FTotal.z * objects[o]._m;

			objects[o].x1.x += deltat * objects[o].FTotal.x * objects[o]._m;
			objects[o].x1.y += deltat * objects[o].FTotal.y * objects[o]._m;
			objects[o].x1.z += deltat * objects[o].FTotal.z * objects[o]._m;

      objects[o].w1.x += deltat * objects[o].MTotal.x * objects[o]._I;
      objects[o].w1.y += deltat * objects[o].MTotal.y * objects[o]._I;
      objects[o].w1.z += deltat * objects[o].MTotal.z * objects[o]._I;



  }


	}

}

// functions for periodic boundary conditions

void DEM::xPbcs(){

	for (int n ; n < particles.size(); n++){


		if (particles[n].x0.x > demSize[0]){

			particles[n].x0.x -= demSize[0];

		}

		if (particles[n].x0.x < 0.0){

			particles[n].x0.x += demSize[0];

		}


	}
}

void DEM::yPbcs(){

	for (int n ; n < particles.size(); n++){


		if (particles[n].x0.y > demSize[1]){

			particles[n].x0.y -= demSize[1];

		}

		if (particles[n].x0.y < 0.0){

			particles[n].x0.y += demSize[1];

		}


	}



}

void DEM::zPbcs(){


	for (int n ; n < particles.size(); n++){


		if (particles[n].x0.z > demSize[2]){

			particles[n].x0.z -= demSize[2];

		}

		if (particles[n].x0.z < 0.0){

			particles[n].x0.z += demSize[2];

		}


	}


}
