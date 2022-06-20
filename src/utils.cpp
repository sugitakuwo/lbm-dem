

#include <vector>
#include <iostream>
//
#include "utils.h"

using namespace std;

/*/////////////////////////////////////////////////////////////////////////////////////////
// WALL
/////////////////////////////////////////////////////////////////////////////////////////*/

double wall::dist(tVect pt) const{ // ACHTUNG!! This assume vector n is unitary
    const tVect d=pt-p;
    return n.dot(d);
}

void wall::wallShow() const {
    cout<<"Wall number "<<index<<": ";
    cout<<"with base point:";
    p.show();
    cout<<" and normal vector:";
    n.show();
    cout<<" moving:"<<moving<<"\n";
    if (moving) {
        cout<<"with translational velocity:";
        vel.show();
        cout<<"\n";
        cout<<"and rotational velocity:";
        omega.show();
        cout<<" around:";
        rotCenter.show();
        cout<<"\n";
    }
}

tVect wall::getSpeed(tVect pt) const {
    if (moving) {
        // distance rotation center
        const tVect distFromCenter=pt-rotCenter;
        // distance from axes
        const tVect distFromAxes=distFromCenter-(distFromCenter.dot(n))*n;
        // tangential velocity
        return vel+omega.cross(distFromAxes);
    }
    else {
        return tVect(0.0,0.0,0.0);
    }
}



/*/////////////////////////////////////////////////////////////////////////////////////////
// TOPOGRAPHY
/////////////////////////////////////////////////////////////////////////////////////////*/

void topography::readFromFile(string& topographyFileName) {
    ifstream topographyFileID;
    cout << "Reading " << topographyFileName.c_str() << endl;
    // open topography file
    topographyFileID.open(topographyFileName.c_str(), ios::in);
    ASSERT(topographyFileID.is_open());
    topographyFileID>>sizeX;
    topographyFileID>>sizeY;
    for (int iy = 0; iy < sizeY; iy++) {
        vecList pointRow;
        doubleList fluidLevelRow;
        for (int ix = 0; ix < sizeX; ix++) {
            double x, y, z, h;
            topographyFileID>>x;
            topographyFileID>>y;
            topographyFileID>>z;
            topographyFileID>>h;
            // add to row
            tVect coordHere(x, y, z);
            pointRow.push_back(coordHere);
            fluidLevelRow.push_back(h);
            // save coordinates
            if (iy == 0) {
                coordX.push_back(x);
            }
            if (ix == 0) {
                coordY.push_back(y);
            }
        }
        points.push_back(pointRow);
        fluidLevel.push_back(fluidLevelRow);
    }
    corner1X=coordX[0];
    corner1Y=coordY[0];
    corner2X=coordX[sizeX-1];
    corner2Y=coordY[sizeY-1];
    spacingX=coordX[1]-coordX[0];
    spacingY=coordY[1]-coordY[0];

}

void topography::show() {

    cout<<"Topography with size "<<sizeX<<" x "<<sizeY<<endl;
    cout<<"Points with spacing "<<spacingX<<" x "<<spacingY<<":"<<endl;
//    for (int iy = 0; iy < sizeY; iy++) {
//        for (int ix = 0; ix < sizeX; ix++) {
//            cout<<"("<<ix<<","<<iy<<")"<<" "<<points[iy][ix].dot(Xp)<<" "<<points[iy][ix].dot(Yp)<<" "<<points[iy][ix].dot(Zp)<<" "<<fluidLevel[iy][ix]<<endl;
//        }
//    }

}


void topography::getReferenceTriangle(const tVect& point, tVect& point1, tVect& point2, tVect& point3, tVect& planeNormal)  const{
    
        // coordinates of the point
    const double pointX = point.dot(Xp);
    const double pointY = point.dot(Yp);

    // approximate indices of the point, projected on the xy plane
    const double indexX = (pointX - corner1X) / spacingX;
    const double indexY = (pointY - corner1Y) / spacingY;

    // indices of the inscribing rectangle
    const int negX = floor(indexX);
    int posX = ceil(indexX);
    if (negX==posX) posX++;
    const int negY = floor(indexY);
    int posY = ceil(indexY);
    if (negY==posY) posY++;
    //cout<<negX<<" "<<posX<<" "<<negY<<" "<<posY<<endl;
    // use the approximate indices to decide which triangle to use
    bool useTop = false;
    if ((indexX - double(negX))<(indexY - double(negY))) {
        useTop = true;
    }

    // setup triangle to use
    point1 = points[negY][negX];
    if (useTop) {
        point2 = points[posY][posX];
        point3 = points[posY][negX];
    } else {
        point2 = points[negY][posX];
        point3 = points[posY][posX];
    }

    // compute triangle normal
    const tVect side1=point2-point1;
    const tVect side2=point3-point1;
    
    planeNormal=side1.cross(side2);

}

double topography::distance(const tVect& point)  const{

    tVect point1,point2,point3;
    tVect planeNormal;
    
    getReferenceTriangle(point,point1,point2,point3,planeNormal);

    // vector connecting the interrogation point and one point of the triangle
    const tVect connection=point-point1;
    // distance is the dot product of these two
    return connection.dot(planeNormal);
    
}

double topography::directionalDistance(const tVect& point, const tVect& dir)  const{
    // distance point to topography, measured along a line defined by dir (normalized vector)
    // see https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
    
    tVect point1,point2,point3;
    tVect planeNormal;
    
    getReferenceTriangle(point,point1,point2,point3,planeNormal);

    // vector connecting the interrogation point and one point of the triangle
    const tVect connection=point-point1;
    // distance is the dot product of these two
    return -1.0*connection.dot(planeNormal)/dir.dot(planeNormal);
    
}


/*/////////////////////////////////////////////////////////////////////////////////////////
// PERIODIC BOUNDARY CONDITION (PBC)
/////////////////////////////////////////////////////////////////////////////////////////*/

void pbc::setPlanes(){
    pl1.p=p;
    pl2.p=p+v;
    pl1.n=v/v.norm();
    pl2.n=-1.0*pl1.n;
}

void pbc::pbcShow() const {
    cout<<"Periodic condition number "<<index<<"\n";
    cout<<"with base point: ";
    p.show();
    cout<<" and translation vector:";
    v.show();
    cout<<"\n";

    cout<<"plane 1: ";
    pl1.p.show();
    cout<<" and normal vector:";
    pl1.n.show();
    cout<<"\n";

    cout<<"plane 2: ";
    pl2.p.show();
    cout<<" and normal vector:";
    pl2.n.show();
    cout<<"\n";

}


