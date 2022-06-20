/*
 * File:   utils.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:11 PM
 */

#ifndef UTILS_H
#define	UTILS_H

#include "myvector.h"
#include "macros.h"
class elmt;

/*//////////////////////////////////////////////////////////////////////////////////////////////////////
// CLASS  DEFINITIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////*/

class wall {
public:
    // index
    unsigned int index;
	//global index to avoid confusion with particles
	unsigned int wallIndex;
    // normal vector
    tVect n;
    // basic point
    tVect p;
    // center of rotation
    tVect rotCenter;
    // rotational speed
    tVect omega;
    // translational speed
    tVect vel;
    // hydraulic force on wall
    tVect FHydro;
    // contact force on wall
    tVect FContact;
    // lubrication force on wall
    tVect FLub;
    // total force on wall
    tVect FTotal;
	// double layer force on wall
	tVect FRep;
    // is it a moving wall? (true = moving; false = fixed);
    bool moving;
    // is it a slipping wall?
    bool slip;
    // is it translating?
    bool translating;
    // translation vector
    tVect trans;
    // limits for walls
    bool limited;
    double xMin,xMax,yMin,yMax,zMin,zMax;
    // default constructor
    wall(){
        n=tVect(1.0,0.0,0.0);
        p=tVect(0.0,0.0,0.0);
        rotCenter=tVect(0.0,0.0,0.0);
        omega=tVect(0.0,0.0,0.0);
        vel=tVect(0.0,0.0,0.0);
        FHydro=tVect(0.0,0.0,0.0);
        FContact=tVect(0.0,0.0,0.0);
        FLub=tVect(0.0,0.0,0.0);
        FTotal=tVect(0.0,0.0,0.0);
		FRep=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
        translating=false;
        trans=tVect(0.0,0.0,0.0);
        limited=false;
        xMin=yMin=zMin=0.0;
        xMax=yMax=zMax=HUGE_VAL;
        index=0;
		wallIndex = 0;
    }
    // better constructor
    wall(tVect ip, tVect in){
        n=in;
        p=ip;
        rotCenter=tVect(0.0,0.0,0.0);
        omega=tVect(0.0,0.0,0.0);
        vel=tVect(0.0,0.0,0.0);
        moving=false;
        slip=false;
        translating=false;
        trans=tVect(0.0,0.0,0.0);
        limited=false;
        xMin=yMin=zMin=0.0;
        xMax=yMax=zMax=HUGE_VAL;
        index=0;
		wallIndex = 0;
    }
    // distance point-wall
    double dist(tVect x) const;
    // show wall characteristics
    void wallShow() const;
    // computes the speed of a point of the wall
    tVect getSpeed(tVect pt) const;
};

class topography {
public:
    // index
    unsigned int sizeX,sizeY;
    // coordinate vectors
    doubleList coordX;
    doubleList coordY;
    vecSet points;
    doubleSet fluidLevel;
    double spacingX,spacingY;
    double corner1X,corner1Y;
    double corner2X,corner2Y;
    topography(){
        coordX.clear();
        coordY.clear();
        points.clear();
        fluidLevel.clear();
    }
    void readFromFile(string& fileName);
    void show();
    // computes the distance sfrom the topographical surface
    void getReferenceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const;
    double distance(const tVect& point) const;
    double directionalDistance(const tVect& point, const tVect& dir) const;
};



class pbc {
 // periodic boundary class, defined by a point and a distance vector between the two planes
public:
    // index
    int index;
    // point of a plane
    tVect p;
    // distance between two planes
    tVect v;
    // two periodic planes (to be defined)
    wall pl1, pl2;
    void pbcShow() const;
    //function for the definition of two planes
    void setPlanes();

    pbc(){
        v=tVect(1.0,0.0,0.0);
        p=tVect(0.0,0.0,0.0);
    }
};




#endif	/* UTILS_H */
