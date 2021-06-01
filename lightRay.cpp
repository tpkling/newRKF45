//
// Developed by Joey Garuti, Tyler Martell, Amanda Coughlin, Anne Marie Pooler & Adam St. Amand
// Based on tpk code for PHYS 422
// Modified to work with lens constructor from PHYS 422 by tpk on April 19, 2017
//
// Reworked by tpk and Michael Peterson in spring 2021, using scaling from nearbylens in github
// adding geodesic deviation vectors
//

#include "lightRay.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

static double b11 = 0;
static double b12 = 0;
static double b13 = 0;
static double b14 = 0;
static double b15 = 0;

static double b21 = 1.0/5.0;
static double b22 = 0.0;
static double b23 = 0.0;
static double b24 = 0.0;
static double b25 = 0.0;

static double b31 = 3.0/40.0;
static double b32 = 9.0/40.0;
static double b33 = 0.0;
static double b34 = 0.0;
static double b35 = 0.0;

static double b41 = 3.0/10.0;
static double b42 = -9.0/10.0;
static double b43 = 6.0/5.0;
static double b44 = 0.0;
static double b45 = 0.0;

static double b51 = -11.0/54.0;
static double b52 = 5.0/2.0;
static double b53 = -70.0/27.0;
static double b54 = 35.0/27.0;
static double b55 = 0.0;

static double b61 = 1631.0/55296.0;
static double b62 = 175.0/512.0;
static double b63 = 575.0/13824.0;
static double b64 = 44275.0/110592.0;
static double b65 = 253.0/4096.0;

static double c1 = 37.0/378.0;
static double c2 = 0.0;
static double c3 = 250.0/621.0;
static double c4 = 125.0/594.0;
static double c5 = 0.0;
static double c6 = 512.0/1771.0;

static double c1s = 2825.0/27648.0;
static double c2s = 0.0;
static double c3s = 18575.0/48384.0;
static double c4s = 13525.0/55296.0;
static double c5s = 277.0/14336.0;
static double c6s = 1.0/4.0;


lightRay::lightRay() {
    // default constructor with starting values
   
//    double c_light = 299792458; // in  m/s
//    double t_s = 4.24989e+17; // in s
//    t_s = t_s*c_light; // in m
    double theta = 15.0*(4.84814e-6); // the "35" is a basic number in arc sec.
    
    av_  = 1.0/(1.0+0.45); // lens at z = 0.45
    apv_ = 0.0; // not including expansion in potential
    coords.t_   = 1.0; // today in scaled time = 1
    coords.x_   = 0.0;
    coords.y_   = 0.0;
    coords.z_   = -0.4171847080; // from github nearbylens - scaled
    coords.vt_  = -1.0; // backwards time directed, not sure if this matters
    coords.vy_  = 0.0;
    coords.vx_  = sin(theta);
    coords.vz_  = cos(theta); // makes the ray go mostly down z axis with a little bend in the x dir.
    coords.Xt_  = 0.0;
    coords.Xx_  = 0.0;
    coords.Xy_  = 0.0;
    coords.Xz_  = 0.0;
    coords.VXt_ = 0.0;
    coords.VXy_ = 0.18786816591e-3; //close, not correct yet. Somthing weird is going on... - Peterson
    coords.VXx_ = 0.5418211e-4; // much smaller than theta, probably good
    coords.VXz_ = 0.0;
}
lightRay::lightRay(double t, double x, double y, double z, double vt, double vx, double vy, double vz) {

    // not currently using
    // set the initial values from whatever is passed in
    coords.t_ = t;
    coords.x_ = x;
    coords.y_ = y;
    coords.z_ = z;
    coords.vt_ = vt;
    coords.vx_ = vx;
    coords.vy_ = vy;
    coords.vz_ = vz;

}

/*  ******************************************************** */


// now define the ODE functions - these version were for including some cosmology, but in code, set to constants

double lightRay::computeVtdot() { //computes vt dot

    return (-2.0*vtv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_))
           /(1.0 + 2.0*phiv_);
}

double lightRay::computeVxdot() { //computes vx dot

    return (+ 2.0*vxv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             - phixv_*(vtv_*vtv_ + vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_))/(1.0-2.0*phiv_); 
}

double lightRay::computeVydot() { //computes vy dot

    return (+ 2.0*vyv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             - phiyv_*(vtv_*vtv_ + vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_))/(1.0-2.0*phiv_); 
}

double lightRay::computeVzdot() { //computes vz dot

    return (+ 2.0*vzv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             - phizv_*(vtv_*vtv_ + vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_))/(1.0-2.0*phiv_); 
}

// also add function for geodesic deviation equation
// Fixed per discussion on Mathematica code notation - Peterson

double lightRay::computeVXtdot(){ // computes the acc in the t component of the geodesic deviation vector
    
    return R0301_*(vzv_*vtv_*Xxv_) + R0302_*(vzv_*vtv_*Xyv_) + R0303_*(-vzv_*vzv_*Xtv_ + vzv_*vtv_*Xzv_) 
            + R0331_*(vzv_*vzv_*Xxv_) + R0332_*(vzv_*vzv_*Xyv_); 
}

double lightRay::computeVXxdot(){ // computes the acc in the x component of the geodesic deviation vector
    
    return R1001_*(vtv_*vtv_*Xxv_) + R1331_*(vzv_*vzv_*Xxv_) + R1031_*(vtv_*vzv_*Xxv_)
            + R1002_*(vtv_*vtv_*Xyv_) + R1030_*(vtv_*vzv_*Xtv_ - vtv_*vtv_*Xzv_) + R1301_*(vzv_*vtv_*Xxv_)
            + R1330_*(vzv_*vzv_*Xtv_ - vzv_*vtv_*Xzv_) + R1332_*(vzv_*vzv_*Xyv_);
}

double lightRay::computeVXydot(){ // computes the acc in the y component of the geodesic deviation vector
    
    return R2002_*(vtv_*vtv_*Xyv_) + R2332_*(vzv_*vzv_*Xyv_) + R2032_*(vtv_*vzv_*Xyv_)
            + R2001_*(vtv_*vtv_*Xxv_) + R2030_*(vtv_*vzv_*Xtv_ - vtv_*vtv_*Xzv_) + R2302_*(vzv_*vtv_*Xyv_)
            + R2330_*(vzv_*vzv_*Xtv_ - vzv_*vtv_*Xzv_) + R2331_*(vzv_*vzv_*Xxv_); 
}

double lightRay::computeVXzdot(){ // computes the acc in the z component of the geodesic deviation vector
    
    return R3030_*(vtv_*vzv_*Xtv_ - vtv_*vtv_*Xzv_) + R3032_*(vtv_*vzv_*Xyv_) + R3031_*(vtv_*vzv_*Xxv_)
            + R3002_*(vtv_*vtv_*Xyv_) + R3001_*(vtv_*vtv_*Xxv_);
}


/*  ********************************************************  */

/*  This next block of code implements Cash-Karp RKF4-5 based on NRC */

void lightRay::setRKF45values(Coords coords, double h, BaseLens* Lens, int part) {

//    cout << part << "  ";
    
    if(part == 0){  // used in computing k1
        tv_  = coords.t_;
        xv_  = coords.x_;
        yv_  = coords.y_;
        zv_  = coords.z_;
        vtv_ = coords.vt_;
        vxv_ = coords.vx_;
        vyv_ = coords.vy_;
        vzv_ = coords.vz_;
        
        Xtv_  = coords.Xt_;
        Xxv_  = coords.Xx_;
        Xyv_  = coords.Xy_;
        Xzv_  = coords.Xz_;
        VXtv_ = coords.VXt_;
        VXxv_ = coords.VXx_;
        VXyv_ = coords.VXy_;
        VXzv_ = coords.VXz_;
//        cout <<" in part 0 " << vxv_ << endl;
        
    } // closes if for RK initial step

    if(part == 1){  // used in computing k2
        tv_  = coords.t_ + b21*k1t_;
        xv_  = coords.x_  + b21*k1x_;
        yv_  = coords.y_  + b21*k1y_;
        zv_  = coords.z_  + b21*k1z_;
        vtv_ = coords.vt_ + b21*k1vt_;
        vxv_ = coords.vx_ + b21*k1vx_;
        vyv_ = coords.vy_ + b21*k1vy_;
        vzv_ = coords.vz_ + b21*k1vz_;
        
        Xtv_  = coords.Xt_  + b21*k1Xt_;
        Xxv_  = coords.Xx_  + b21*k1Xx_;
        Xyv_  = coords.Xy_  + b21*k1Xy_;
        Xzv_  = coords.Xz_  + b21*k1Xz_;
        VXtv_ = coords.VXt_ + b21*k1VXt_;
        VXxv_ = coords.VXx_ + b21*k1VXx_;
        VXyv_ = coords.VXy_ + b21*k1VXy_;
        VXzv_ = coords.VXz_ + b21*k1VXz_;
//        cout <<" in part 1 " << vxv_ << endl;
        
    } // closes if for RK 1st step

    if(part == 2){  // used in computing k3
        tv_  = coords.t_  + b41*k1t_ + b42*k2t_ + b43*k3t_;
        xv_  = coords.x_  + b41*k1x_ + b42*k2x_ + b43*k3x_;
        yv_  = coords.y_  + b41*k1y_ + b42*k2y_ + b43*k3y_;
        zv_  = coords.z_  + b41*k1z_ + b42*k2z_ + b43*k3z_;
        vtv_ = coords.vt_ + b41*k1vt_ + b42*k2vt_ + b43*k3vt_;
        vxv_ = coords.vx_ + b41*k1vx_ + b42*k2vx_ + b43*k3vx_;
        vyv_ = coords.vy_ + b41*k1vy_ + b42*k2vy_ + b43*k3vy_;
        vzv_ = coords.vz_ + b41*k1vz_ + b42*k2vz_ + b43*k3vz_;
        
        Xtv_  = coords.Xt_  + b41*k1Xt_ + b42*k2Xt_ + b43*k3Xt_;
        Xxv_  = coords.Xx_  + b41*k1Xx_ + b42*k2Xx_ + b43*k3Xx_;
        Xyv_  = coords.Xy_  + b41*k1Xy_ + b42*k2Xy_ + b43*k3Xy_;
        Xzv_  = coords.Xz_  + b41*k1Xz_ + b42*k2Xz_ + b43*k3Xz_;
        VXtv_ = coords.VXt_ + b41*k1VXt_ + b42*k2VXt_ + b43*k3VXt_;
        VXxv_ = coords.VXx_ + b41*k1VXx_ + b42*k2VXx_ + b43*k3VXx_;
        VXyv_ = coords.VXy_ + b41*k1VXy_ + b42*k2VXy_ + b43*k3VXy_;
        VXzv_ = coords.VXz_ + b41*k1VXz_ + b42*k2VXz_ + b43*k3VXz_;
//        cout <<" in part 2 " << vxv_ << endl;
        
    } // closes if for RK  2nd step

    if(part == 3){  // used in computing k4
        tv_  = coords.t_  + b41*k1t_ + b42*k2t_ + b43*k3t_;
        xv_  = coords.x_  + b41*k1x_ + b42*k2x_ + b43*k3x_;
        yv_  = coords.y_  + b41*k1y_ + b42*k2y_ + b43*k3y_;
        zv_  = coords.z_  + b41*k1z_ + b42*k2z_ + b43*k3z_;
        vtv_ = coords.vt_ + b41*k1vt_ + b42*k2vt_ + b43*k3vt_;
        vxv_ = coords.vx_ + b41*k1vx_ + b42*k2vx_ + b43*k3vx_;
        vyv_ = coords.vy_ + b41*k1vy_ + b42*k2vy_ + b43*k3vy_;
        vzv_ = coords.vz_ + b41*k1vz_ + b42*k2vz_ + b43*k3vz_;
        
        Xtv_  = coords.Xt_  + b41*k1Xt_ + b42*k2Xt_ + b43*k3Xt_;
        Xxv_  = coords.Xx_  + b41*k1Xx_ + b42*k2Xx_ + b43*k3Xx_;
        Xyv_  = coords.Xy_  + b41*k1Xy_ + b42*k2Xy_ + b43*k3Xy_;
        Xzv_  = coords.Xz_  + b41*k1Xz_ + b42*k2Xz_ + b43*k3Xz_;
        VXtv_ = coords.VXt_ + b41*k1VXt_ + b42*k2VXt_ + b43*k3VXt_;
        VXxv_ = coords.VXx_ + b41*k1VXx_ + b42*k2VXx_ + b43*k3VXx_;
        VXyv_ = coords.VXy_ + b41*k1VXy_ + b42*k2VXy_ + b43*k3VXy_;
        VXzv_ = coords.VXz_ + b41*k1VXz_ + b42*k2VXz_ + b43*k3VXz_;
//        cout <<" in part 3 " << vxv_ << "  " << k1vx_ <<"  "<< k2vx_ <<"  " << k3vx_ << endl;
        
    } // closes if for RK  3rd step

    if(part == 4){  // used in computing k5
        tv_  = coords.t_  + b51*k1t_ + b52*k2t_ + b53*k3t_ + b54*k4t_;
        xv_  = coords.x_  + b51*k1x_ + b52*k2x_ + b53*k3x_ + b54*k4x_;
        yv_  = coords.y_  + b51*k1y_ + b52*k2y_ + b53*k3y_ + b54*k4y_;
        zv_  = coords.z_  + b51*k1z_ + b52*k2z_ + b53*k3z_ + b54*k4z_;
        vtv_ = coords.vt_ + b51*k1vt_ + b52*k2vt_ + b53*k3vt_ + b54*k4vt_;
        vxv_ = coords.vx_ + b51*k1vx_ + b52*k2vx_ + b53*k3vx_ + b54*k4vx_;
        vyv_ = coords.vy_ + b51*k1vy_ + b52*k2vy_ + b53*k3vy_ + b54*k4vy_;
        vzv_ = coords.vz_ + b51*k1vz_ + b52*k2vz_ + b53*k3vz_ + b54*k4vz_;
        
        Xtv_  = coords.Xt_  + b51*k1Xt_ + b52*k2Xt_ + b53*k3Xt_ + b54*k4Xt_;
        Xxv_  = coords.Xx_  + b51*k1Xx_ + b52*k2Xx_ + b53*k3Xx_ + b54*k4Xx_;
        Xyv_  = coords.Xy_  + b51*k1Xy_ + b52*k2Xy_ + b53*k3Xy_ + b54*k4Xy_;
        Xzv_  = coords.Xz_  + b51*k1Xz_ + b52*k2Xz_ + b53*k3Xz_ + b54*k4Xz_;
        VXtv_ = coords.VXt_ + b51*k1VXt_ + b52*k2VXt_ + b53*k3VXt_ + b54*k4VXt_;
        VXxv_ = coords.VXx_ + b51*k1VXx_ + b52*k2VXx_ + b53*k3VXx_ + b54*k4VXx_;
        VXyv_ = coords.VXy_ + b51*k1VXy_ + b52*k2VXy_ + b53*k3VXy_ + b54*k4VXy_;
        VXzv_ = coords.VXz_ + b51*k1VXz_ + b52*k2VXz_ + b53*k3VXz_ + b54*k4VXz_;
//        cout <<" in part 4 " << vxv_ << endl;
        
    } // closes if for RK  4th step

    if(part == 5){  // used in computing k6
        tv_  = coords.t_  + b61*k1t_ + b62*k2t_ + b63*k3t_ + b64*k4t_ + b65*k5t_;
        xv_  = coords.x_  + b61*k1x_ + b62*k2x_ + b63*k3x_ + b64*k4x_ + b65*k5x_;
        yv_  = coords.y_  + b61*k1y_ + b62*k2y_ + b63*k3y_ + b64*k4y_ + b65*k5y_;
        zv_  = coords.z_  + b61*k1z_ + b62*k2z_ + b63*k3z_ + b64*k4z_ + b65*k5z_;
        vtv_ = coords.vt_ + b61*k1vt_ + b62*k2vt_ + b63*k3vt_ + b64*k4vt_ + b65*k5vt_;
        vxv_ = coords.vx_ + b61*k1vx_ + b62*k2vx_ + b63*k3vx_ + b64*k4vx_ + b65*k5vx_;
        vyv_ = coords.vy_ + b61*k1vy_ + b62*k2vy_ + b63*k3vy_ + b64*k4vy_ + b65*k5vy_;
        vzv_ = coords.vz_ + b61*k1vz_ + b62*k2vz_ + b63*k3vz_ + b64*k4vz_ + b65*k5vz_;
        
        Xtv_  = coords.Xt_  + b61*k1Xt_ + b62*k2Xt_ + b63*k3Xt_ + b64*k4Xt_ + b65*k5Xt_;
        Xxv_  = coords.Xx_  + b61*k1Xx_ + b62*k2Xx_ + b63*k3Xx_ + b64*k4Xx_ + b65*k5Xx_;
        Xyv_  = coords.Xy_  + b61*k1Xy_ + b62*k2Xy_ + b63*k3Xy_ + b64*k4Xy_ + b65*k5Xy_;
        Xzv_  = coords.Xz_  + b61*k1Xz_ + b62*k2Xz_ + b63*k3Xz_ + b64*k4Xz_ + b65*k5Xz_;
        VXtv_ = coords.VXt_ + b61*k1VXt_ + b62*k2VXt_ + b63*k3VXt_ + b64*k4VXt_ + b65*k5VXt_;
        VXxv_ = coords.VXx_ + b61*k1VXx_ + b62*k2VXx_ + b63*k3VXx_ + b64*k4VXx_ + b65*k5VXx_;
        VXyv_ = coords.VXy_ + b61*k1VXy_ + b62*k2VXy_ + b63*k3VXy_ + b64*k4VXy_ + b65*k5VXy_;
        VXzv_ = coords.VXz_ + b61*k1VXz_ + b62*k2VXz_ + b63*k3VXz_ + b64*k4VXz_ + b65*k5VXz_;
//        cout <<" in part 5 " << vxv_ << endl;
        
    } // closes if for RK  5th step

//    av_ = 1.0;   // if not static, update these here
//    apv_ = 0.0;  // for static values, set in constructor
//    rv_ = sqrt(xv_*xv_ + yv_*yv_ + zv_*zv_); // this being used?

    phiv_   = Lens->getU(xv_, yv_, zv_);
    phixv_  = Lens->getUdx(xv_, yv_, zv_, h); // passing h here, in the numerical derivative, defines d = h/10 and moves that far
    phiyv_  = Lens->getUdy(xv_, yv_, zv_, h);
    phizv_  = Lens->getUdz(xv_, yv_, zv_, h);
    phixxv_ = Lens->getUdxx(xv_, yv_, zv_, h); // set second deriv values
    phiyyv_ = Lens->getUdyy(xv_, yv_, zv_, h);
    phizzv_ = Lens->getUdzz(xv_, yv_, zv_, h);
    phixyv_ = Lens->getUdxy(xv_, yv_, zv_, h);
    phixzv_ = Lens->getUdxz(xv_, yv_, zv_, h);
    phiyzv_ = Lens->getUdyz(xv_, yv_, zv_, h);

    R0301_ = - phixzv_; 
    R0302_ = - phiyzv_;
    R0303_ = - phizzv_; // - phittv_ working with no time dependence at the moment
    R0331_ = 0; // + phitxv_
    R0332_ = 0; // + phityv_
    R1001_ = - phixxv_ ; // - phittv_
    R1331_ = - phixxv_ - phizzv_;
    R1031_ = 0; // - phitzv_
    R1002_ = - phixyv_;
    R1030_ = - phixzv_;
    R1301_ = 0; // - phitzv_
    R1330_ = 0; // - phitxv_
    R1332_ = - phixyv_;
    R2002_ = - phiyyv_; // - phittv_
    R2332_ = - phiyyv_ - phizzv_;
    R2032_ = 0; // - phitzv_
    R2001_ = - phixyv_;
    R2030_ = 0; // - phitzv_
    R2302_ = 0; // - phitzv_
    R2330_ = 0; // - phityv_
    R2331_ = - phixyv_;
    R3030_ = - phizzv_; // - phittv_
    R3032_ = 0; // - phityv_
    R3031_ = 0; // - phitxv_
    R3002_ = 0; // - phiyzv_
    R3001_ = - phixzv_;
    
} // closes function to set values

/* This file takes 1 RKF step */

//void lightRay::takeRKF45step(Coords coords, Coords &newCoords, double errvec[], double h, BaseLens* Lens) {
Coords lightRay::takeRKF45step(Coords coords, double errvec[], double h, BaseLens* Lens){    
    setRKF45values(coords, h, Lens, 0); // setting up initial values

//    cout << "0: " << vxv_ << "  " << vzv_ << endl;
    
    k1t_ = h*vtv_;
    k1x_ = h*vxv_;
    k1y_ = h*vyv_;
    k1z_ = h*vzv_;
    k1vt_ = h*computeVtdot();
    k1vx_ = h*computeVxdot();
    k1vy_ = h*computeVydot();
    k1vz_ = h*computeVzdot();
    
    k1Xt_ = h*VXtv_;
    k1Xx_ = h*VXxv_;
    k1Xy_ = h*VXyv_;
    k1Xz_ = h*VXzv_;
    k1VXt_ = h*computeVXtdot();
    k1VXx_ = h*computeVXxdot();
    k1VXy_ = h*computeVXydot();
    k1VXz_ = h*computeVXzdot();

//    setRKF45values(coords, h, Lens, 1); // set up values to use in k2's

//    cout << "1: " << vxv_ << "  " << vzv_ << endl;
    k2t_ = h*vtv_;
    k2x_ = h*vxv_;
    k2y_ = h*vyv_;
    k2z_ = h*vzv_;
    k2vt_ = h*computeVtdot();
    k2vx_ = h*computeVxdot();
    k2vy_ = h*computeVydot();
    k2vz_ = h*computeVzdot();

    k2Xt_ = h*VXtv_;
    k2Xx_ = h*VXxv_;
    k2Xy_ = h*VXyv_;
    k2Xz_ = h*VXzv_;
    k2VXt_ = h*computeVXtdot();
    k2VXx_ = h*computeVXxdot();
    k2VXy_ = h*computeVXydot();
    k2VXz_ = h*computeVXzdot();

//    setRKF45values(coords, h, Lens, 2); // set up values to use in k3's

//    cout << "2: " << vxv_ << "  " << vzv_ << endl;
    k3t_ = h*vtv_;
    k3x_ = h*vxv_;
    k3y_ = h*vyv_;
    k3z_ = h*vzv_;
    k3vt_ = h*computeVtdot();
    k3vx_ = h*computeVxdot();
    k3vy_ = h*computeVydot();
    k3vz_ = h*computeVzdot();
    
    k3Xt_ = h*VXtv_;
    k3Xx_ = h*VXxv_;
    k3Xy_ = h*VXyv_;
    k3Xz_ = h*VXzv_;
    k3VXt_ = h*computeVXtdot();
    k3VXx_ = h*computeVXxdot();
    k3VXy_ = h*computeVXydot();
    k3VXz_ = h*computeVXzdot();

//    setRKF45values(coords, h, Lens, 3); // set up values to use in k4's

//    cout << "3: " << vxv_ << "  " << vzv_ << endl;
    k4t_ = h*vtv_;
    k4x_ = h*vxv_;
    k4y_ = h*vyv_;
    k4z_ = h*vzv_;
    k4vt_ = h*computeVtdot();
    k4vx_ = h*computeVxdot();
    k4vy_ = h*computeVydot();
    k4vz_ = h*computeVzdot();

    k4Xt_ = h*VXtv_;
    k4Xx_ = h*VXxv_;
    k4Xy_ = h*VXyv_;
    k4Xz_ = h*VXzv_;
    k4VXt_ = h*computeVXtdot();
    k4VXx_ = h*computeVXxdot();
    k4VXy_ = h*computeVXydot();
    k4VXz_ = h*computeVXzdot();

    setRKF45values(coords, h, Lens, 4); // set up values to use in k5's

//    cout << "4: " << vxv_ << "  " << vzv_ << endl;
    k5t_ = h*vtv_;
    k5x_ = h*vxv_;
    k5y_ = h*vyv_;
    k5z_ = h*vzv_;
    k5vt_ = h*computeVtdot();
    k5vx_ = h*computeVxdot();
    k5vy_ = h*computeVydot();
    k5vz_ = h*computeVzdot();
    
    k5Xt_ = h*VXtv_;
    k5Xx_ = h*VXxv_;
    k5Xy_ = h*VXyv_;
    k5Xz_ = h*VXzv_;
    k5VXt_ = h*computeVXtdot();
    k5VXx_ = h*computeVXxdot();
    k5VXy_ = h*computeVXydot();
    k5VXz_ = h*computeVXzdot();

    setRKF45values(coords, h, Lens, 5); // set up values to use in k6's

//    cout << "5: " << vxv_ << "  " << vzv_ << endl;
    k6t_ = h*vtv_;
    k6x_ = h*vxv_;
    k6y_ = h*vyv_;
    k6z_ = h*vzv_;
    k6vt_ = h*computeVtdot();
    k6vx_ = h*computeVxdot();
    k6vy_ = h*computeVydot();
    k6vz_ = h*computeVzdot();

    k6Xt_ = h*VXtv_;
    k6Xx_ = h*VXxv_;
    k6Xy_ = h*VXyv_;
    k6Xz_ = h*VXzv_;
    k6VXt_ = h*computeVXtdot();
    k6VXx_ = h*computeVXxdot();
    k6VXy_ = h*computeVXydot();
    k6VXz_ = h*computeVXzdot();

    // save values into temporary spots
    
//   cout << k1x_ <<" " << k2x_ <<" " << k3x_ <<" " << k4x_ <<" " << k5x_ <<" " << k6x_ << endl;
//    char w1;
//    cin >> w1;

    Coords newCoords;
    
    newCoords.t_ = coords.t_ + c1*k1t_ + c2*k2t_ + c3*k3t_ + c4*k4t_ + c5*k5t_ + c6*k6t_;
    newCoords.x_ = coords.x_ + c1*k1x_ + c2*k2x_ + c3*k3x_ + c4*k4x_ + c5*k5x_ + c6*k6x_;
    newCoords.y_ = coords.y_ + c1*k1y_ + c2*k2y_ + c3*k3y_ + c4*k4y_ + c5*k5y_ + c6*k6y_;
    newCoords.z_ = coords.z_ + c1*k1z_ + c2*k2z_ + c3*k3z_ + c4*k4z_ + c5*k5z_ + c6*k6z_;
    newCoords.vt_ = coords.vt_ + c1*k1vt_ + c2*k2vt_ + c3*k3vt_ + c4*k4vt_ + c5*k5vt_ + c6*k6vt_;
    newCoords.vx_ = coords.vx_ + c1*k1vx_ + c2*k2vx_ + c3*k3vx_ + c4*k4vx_ + c5*k5vx_ + c6*k6vx_;
    newCoords.vy_ = coords.vy_ + c1*k1vy_ + c2*k2vy_ + c3*k3vy_ + c4*k4vy_ + c5*k5vy_ + c6*k6vy_;
    newCoords.vz_ = coords.vz_ + c1*k1vz_ + c2*k2vz_ + c3*k3vz_ + c4*k4vz_ + c5*k5vz_ + c6*k6vz_;
    
    newCoords.Xt_ = coords.Xt_ + c1*k1Xt_ + c2*k2Xt_ + c3*k3Xt_ + c4*k4Xt_ + c5*k5Xt_ + c6*k6Xt_;
    newCoords.Xx_ = coords.Xx_ + c1*k1Xx_ + c2*k2Xx_ + c3*k3Xx_ + c4*k4Xx_ + c5*k5Xx_ + c6*k6Xx_;
    newCoords.Xy_ = coords.Xy_ + c1*k1Xy_ + c2*k2Xy_ + c3*k3Xy_ + c4*k4Xy_ + c5*k5Xy_ + c6*k6Xy_;
    newCoords.Xz_ = coords.Xz_ + c1*k1Xz_ + c2*k2Xz_ + c3*k3Xz_ + c4*k4Xz_ + c5*k5Xz_ + c6*k6Xz_;
    newCoords.VXt_ = coords.VXt_ + c1*k1VXt_ + c2*k2VXt_ + c3*k3VXt_ + c4*k4VXt_ + c5*k5VXt_ + c6*k6VXt_;
    newCoords.VXx_ = coords.VXx_ + c1*k1VXx_ + c2*k2VXx_ + c3*k3VXx_ + c4*k4VXx_ + c5*k5VXx_ + c6*k6VXx_;
    newCoords.VXy_ = coords.VXy_ + c1*k1VXy_ + c2*k2VXy_ + c3*k3VXy_ + c4*k4VXy_ + c5*k5VXy_ + c6*k6VXy_;
    newCoords.VXz_ = coords.VXz_ + c1*k1VXz_ + c2*k2VXz_ + c3*k3VXz_ + c4*k4VXz_ + c5*k5VXz_ + c6*k6VXz_;
    
    // store err values

    errvec[0] = fabs((c1-c1s)*k1t_ + (c2-c2s)*k2t_ + (c3-c3s)*k3t_ + (c4-c4s)*k4t_ + (c5-c5s)*k5t_ + (c6-c6s)*k6t_);
    errvec[1] = fabs((c1-c1s)*k1x_ + (c2-c2s)*k2x_ + (c3-c3s)*k3x_ + (c4-c4s)*k4x_ + (c5-c5s)*k5x_ + (c6-c6s)*k6x_);
    errvec[2] = fabs((c1-c1s)*k1y_ + (c2-c2s)*k2y_ + (c3-c3s)*k3y_ + (c4-c4s)*k4y_ + (c5-c5s)*k5y_ + (c6-c6s)*k6y_);
    errvec[3] = fabs((c1-c1s)*k1z_ + (c2-c2s)*k2z_ + (c3-c3s)*k3z_ + (c4-c4s)*k4z_ + (c5-c5s)*k5z_ + (c6-c6s)*k6z_);
    errvec[4] = fabs((c1-c1s)*k1vt_ + (c2-c2s)*k2vt_ + (c3-c3s)*k3vt_ + (c4-c4s)*k4vt_ + (c5-c5s)*k5vt_ + (c6-c6s)*k6vt_);
    errvec[5] = fabs((c1-c1s)*k1vx_ + (c2-c2s)*k2vx_ + (c3-c3s)*k3vx_ + (c4-c4s)*k4vx_ + (c5-c5s)*k5vx_ + (c6-c6s)*k6vx_);
    errvec[6] = fabs((c1-c1s)*k1vy_ + (c2-c2s)*k2vy_ + (c3-c3s)*k3vy_ + (c4-c4s)*k4vy_ + (c5-c5s)*k5vy_ + (c6-c6s)*k6vy_);
    errvec[7] = fabs((c1-c1s)*k1vz_ + (c2-c2s)*k2vz_ + (c3-c3s)*k3vz_ + (c4-c4s)*k4vz_ + (c5-c5s)*k5vz_ + (c6-c6s)*k6vz_);
    errvec[8] = fabs((c1-c1s)*k1Xt_ + (c2-c2s)*k2Xt_ + (c3-c3s)*k3Xt_ + (c4-c4s)*k4Xt_ + (c5-c5s)*k5Xt_ + (c6-c6s)*k6Xt_);
    errvec[9] = fabs((c1-c1s)*k1Xx_ + (c2-c2s)*k2Xx_ + (c3-c3s)*k3Xx_ + (c4-c4s)*k4Xx_ + (c5-c5s)*k5Xx_ + (c6-c6s)*k6Xx_);
    errvec[10] = fabs((c1-c1s)*k1Xy_ + (c2-c2s)*k2Xy_ + (c3-c3s)*k3Xy_ + (c4-c4s)*k4Xy_ + (c5-c5s)*k5Xy_ + (c6-c6s)*k6Xy_);
    errvec[11] = fabs((c1-c1s)*k1Xz_ + (c2-c2s)*k2Xz_ + (c3-c3s)*k3Xz_ + (c4-c4s)*k4Xz_ + (c5-c5s)*k5Xz_ + (c6-c6s)*k6Xz_);
    errvec[12] = fabs((c1-c1s)*k1VXt_ + (c2-c2s)*k2VXt_ + (c3-c3s)*k3VXt_ + (c4-c4s)*k4VXt_ + (c5-c5s)*k5VXt_ + (c6-c6s)*k6VXt_);
    errvec[13] = fabs((c1-c1s)*k1VXx_ + (c2-c2s)*k2VXx_ + (c3-c3s)*k3VXx_ + (c4-c4s)*k4VXx_ + (c5-c5s)*k5VXx_ + (c6-c6s)*k6VXx_);
    errvec[14] = fabs((c1-c1s)*k1VXy_ + (c2-c2s)*k2VXy_ + (c3-c3s)*k3VXy_ + (c4-c4s)*k4VXy_ + (c5-c5s)*k5VXy_ + (c6-c6s)*k6VXy_);
    errvec[15] = fabs((c1-c1s)*k1VXz_ + (c2-c2s)*k2VXz_ + (c3-c3s)*k3VXz_ + (c4-c4s)*k4VXz_ + (c5-c5s)*k5VXz_ + (c6-c6s)*k6VXz_);
    
    return newCoords;

}

void lightRay::runRay(BaseLens* Lens) {
    double h = 1.0e-5;
    double hnew;
    int countSafe = 0;
    double zstop = 0.25; // in scaled units
    double errmax;
    double SAFETY = 0.9;
    double EPS = 1.0e-11;
    double errvec[16];
    double s=0;
    Coords newCoords;
    ofstream fout;
    fout.open("path.csv", ios::app);
    
    
    while (coords.z_ < zstop) {
        countSafe++;
        errmax = 0.0;

        newCoords = takeRKF45step(coords, errvec, h, Lens); // runs a ray for one step
        for (int i = 0; i < 16; i++) {
            if (errvec[i] > errmax) {
                errmax = errvec[i];
            }
        }

        
        if (EPS / errmax > 1.0) {
            coords = newCoords;
            s= s+h;
            hnew = SAFETY*h*pow((EPS/errmax), 0.25);
            if (hnew>1.1*h) {
                h = 1.1*h;
            }
            if (hnew > 1.0e-3)
                h = 1.0e-3;
            else h = hnew;

            //std::cout << "h get bigger = " << h << std::endl;
        }
        else {
            hnew = SAFETY*h*pow(EPS/errmax, 0.2);
            if(hnew < 0.7*h){
                h = 0.7*h;
            }
            else h = hnew;
            //std::cout << "h get smaller = " << h << std::endl;
        }

       if (countSafe % 40 == 0) {
            cout << countSafe << " , " << s <<" , "  << coords.x_ << ",  " << coords.z_ << " , " << 
                coords.vx_ << ",  " << coords.vz_ << " , " <<coords.Xx_ << " , " << coords.Xy_ << " , " << coords.Xz_ << endl;
            fout << countSafe << " , " << s <<" , "  << coords.x_ << ",  " << coords.z_ << " , " << 
                coords.vx_ << ",  " << coords.vz_ << " , " <<coords.Xx_ << " , " << coords.Xy_ << " , " << coords.Xz_ << endl;
        }
        if(countSafe > 20000) {
            cout << "breaking steps" << endl;
            break;
        }
    }//closes while loop
    
    fout.close();
    
}

/* OLD CODE */

/*
double lightRay::computeVtdot() { //computes t dot

    return (-av_*apv_*(1.0-2.0*phiv_)*(vxv_*vxv_+vyv_*vyv_ + vzv_*vzv_)
            -2.0*vtv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_))
           /(1.0 + 2.0*phiv_);
}

double lightRay::computeVxdot() { //computes x dot

    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vxv_
             + 2.0*av_*av_*vxv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phixv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); 
}

double lightRay::computeVydot() { //computes y dot
    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vyv_
             + 2.0*av_*av_*vyv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phiyv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); 
}

double lightRay::computeVzdot() { //computes z dot

    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vzv_
             + 2.0*av_*av_*vzv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phizv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); 
}


*/