/*
 * File:   node.h
 * Author: aleonard
 *
 * Created on November 5, 2013, 4:10 PM
 */

#ifndef NODE_H
#define	NODE_H

#include <stdlib.h>

#include "lattice.h"
#include "myvector.h"


enum Rheology {NEWTONIAN, BINGHAM, FRICTIONAL, VOELLMY};

class FluidMaterial{
public:

    // initial density
    double initDensity;
    // initial viscosity
    double initDynVisc;
    // rheology model
    Rheology rheologyModel;
    // parameters for Bingham fluid (rheologyModel=BINGHAM)
    double yieldStress, plasticVisc;    
    // parameter for friction (rheologyModel=FRICTIONAL,VOELLMY)
    double frictionCoefFluid;
    // constant for collisional/turbulent term (rheologyModel=VOELLMY)
    double voellmyCoefficient;
    // parameter for collisional/turbulent term (rheologyModel=VOELLMY)
    double voellmyConstant;
    // activation of Smagorinsky turbulence model
    bool turbulenceOn;
    // parameter for Smagorinsky turbulence model
    double turbConst;
            // default constructor
    FluidMaterial(){
        initDensity=1.0;
        initDynVisc=(1.0-0.5)/3.0;
        rheologyModel=NEWTONIAN;
        yieldStress=0.0;
        plasticVisc=0.0;
        frictionCoefFluid=0.0;
        voellmyCoefficient=0.0;
        turbulenceOn=false;
        turbConst=0.0;
    }
};

class node{
    // lattice node declaration
public:
    // probability density functions (one per lattice direction)
    double f[lbmDirec];
    // support variables for streaming
    double fs[lbmDirec];
    // density in the node
    double n;
    // velocity of fluid in node
    tVect u;
    tVect hydroForce;
    // mass functions
    double mass, newMass;
    // viscosity
    double visc;
    // neighbor nodes index
    unsigned int d[lbmDirec];
    // informations for curved boundaries
    curve* curved;
	

	
    // default constructor
    node(){
        n=0.0;
        u.reset();
        hydroForce.reset();
        newMass=mass=0.0;
        f[0]=f[1]=f[2]=f[3]=f[4]=f[5]=f[6]=f[7]=f[8]=f[9]=0.0;
        f[10]=f[11]=f[12]=f[13]=f[14]=f[15]=f[16]=f[17]=f[18]=0.0;
        fs[0]=fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=fs[6]=fs[7]=fs[8]=fs[9]=0.0;
        fs[10]=fs[11]=fs[12]=fs[13]=fs[14]=fs[15]=fs[16]=fs[17]=fs[18]=0.0;
        d[0]=d[1]=d[2]=d[3]=d[4]=d[5]=d[6]=d[7]=d[8]=d[9]=0;
        d[10]=d[11]=d[12]=d[13]=d[14]=d[15]=d[16]=d[17]=d[18]=0;
        curved=0;
    }
    // getting liquid fraction
    double liquidFraction() const;
    // shows the population in case debugging requires it
    void showPopulations() const;
    // functions working on distribution functions
    void initialize(const double& density,const tVect& velocity, const double& massFunction, const double& viscosity, const tVect& lbF);
    void restart(const double& restDensity, const tVect& restVelocity, const double& restMassFunction, const double& restViscosity, const double restF[]) ;
    void setEquilibrium(const tVect& F);
    void reconstruct();
    void shiftVelocity(const tVect& F);
    void computeEquilibrium(double feq[]);
    //void computeShearRate(const double feq[], const double& turbConst, const double& plasticVisc, const double& yieldStress);
    void computeApparentViscosity(const double feq[], const FluidMaterial& fluidMaterial);
    void resetViscosity(const double& initViscosity);
    void solveCollision(const double feq[]);
    void addForce(const tVect& F);
    void collideFNN(const double& turbConst, const tVect& F, const double& plasticVisc, const double& yieldStress);
    void collideNN(const double& turbConst, const double& plasticVisc, const double& yieldStress);
    void collideF(const tVect& F, const double& initVisc);
    void collide(const double& initVisc);
    void store();
    tVect bounceBackForce(const unsigned int& j, const double staticPres[], const double& BBi) const;
    double massStream(const unsigned int& sdir) const;
};

class nodeType{
private:
    // identification of node type
    // 0 = liquid cells
    // 1 = solid particles
    // 2 = gas (=void) cells
    // 3 = interface
    // 4 = periodic
    // 5 = slip stationary planar walls
    // 6 = slip moving planar walls
    // 7 = no-slip stationary planar walls
    // 8 = no-slip moving planar walls
    // 9 = cylinder walls
    // 10 = object walls
    // 11 = topographic walls
    unsigned char t;
    // particle flag
    bool p;
    // solid index
    unsigned short int solidIndex;
public:
        // functions for identification of type
    bool isActive() const;
    bool isWall() const;
    bool isCurvedWall() const;
    bool isSlipStatWall() const; 
    bool isSlipDynWall() const;
    bool isStatWall() const;
    bool isDynWall() const;
    bool isCylinderWall() const;
    bool isObjectWall() const;
    bool isTopography() const;
    bool isParticle() const;
    bool isInterface() const;
    bool isFluid() const;
    bool isGas() const;
    bool isPeriodic() const;
    // functions for change of type
    void setParticle();
    void setSlipStatWall();
    void setSlipDynWall();
    void setStatWall();
    void setDynWall();
    void setCylinderWall();
    void setObjectWall();
    void setTopography();
    void setInterface();
    void setFluid();
    void setGas();
    void setPeriodic();
    void setType(unsigned int& typ);
    unsigned int getType() const;
    // particle flag
    bool isInsideParticle() const;
    void setInsideParticle();
    void setOutsideParticle();
    // functions to get and assign index
    unsigned int getSolidIndex() const;
    void setSolidIndex(const unsigned int& ind);
};

class curve{
public:
    double delta[lbmDirec];
    double m1[lbmDirec], m2[lbmDirec];
    curve(){
        for (int j=1; j<lbmDirec; ++j) {
            delta[j]=0.0;
            m1[j]=0.0;
            m2[j]=0.0;
        }
    }
    void computeCoefficients();
    double getChi(const unsigned int& j, const double& tau) const;
};

class measureUnits{
public:
    // conversion units /////////////////////////////////////////////////////////////////////////
    // fluid and granular matter are solved in different measure units
    // for a reference, check Feng, Han, Owen, 2007
    // unit length
    double Length;
    //unit time
    double Time;
    //unit density
    double Density;
    // secondary units: speed, acceleration, mass, force,
    double Volume, Speed, Accel, AngVel, Force, Torque, Mass, KinVisc, DynVisc, Stress, Pressure;
    // inverted unit length for efficiency
    double invLength;
    // constructor
    measureUnits(){
        Volume=1.0;
        Length=1.0;
        Time=1.0;
        Density=1.0;
        Speed=1.0;
        Accel=1.0;
        AngVel=1.0;
        Force=1.0;
        Torque=1.0;
        Mass=1.0;
        KinVisc=1.0;
        DynVisc=1.0;
        Stress=1.0;
        Pressure=1.0;
    }
    void setComposite();
};

#endif	/* NODE_H */

