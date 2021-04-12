//
// Developed by Joey Garuti, Tyler Martell, Amanda Coughlin, Anne Marie Pooler & Adam St. Amand
// Based on tpk code for PHYS 422
// Modified to work with lens constructor from PHYS 422 by tpk on April 19, 2017
//

#include "lightRay.h"
#include <cmath>

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
    // default constructor with dummy values

    t_ = 0.0;
    x_ = 0.0;
    y_ = 0.0;
    z_ = 20.0;
    vt_ = 1.0;
    vy_ = 0.0;
    vx_ = sin(0.3);
    vz_ = cos(0.3); // makes the ray go mostly down z axis with a little bend in the x dir.
}
lightRay::lightRay(double t, double x, double y, double z, double vt, double vx, double vy, double vz) {

    // set the initial values from whatever is passed in
    t_ = t;
    x_ = x;
    y_ = y;
    z_ = z;
    vt_ = vt;
    vx_ = vx;
    vy_ = vy;
    vz_ = vz;

}

/*  ******************************************************** */


// now define the ODE functions


double lightRay::computeVtdot() { //computes t dot

    return (-av_*apv_*(1.0-2.0*phiv_)*(vxv_*vxv_+vyv_*vyv_ + vzv_*vzv_)
            -2.0*vtv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_))
           /(1.0 + 2.0*phiv_);
}

double lightRay::computeVxdot() { //computes x dot

    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vxv_
             + 2.0*av_*av_*vxv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phixv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); // replace with stuff
}

double lightRay::computeVydot() { //computes y dot
    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vyv_
             + 2.0*av_*av_*vyv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phiyv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); // replace with stuff
}

double lightRay::computeVzdot() { //computes z dot

    return -(2.0*av_*apv_*(1-2.0*phiv_)*vtv_*vzv_
             + 2.0*av_*av_*vzv_*(phixv_*vxv_ + phiyv_*vyv_ + phizv_*vzv_)
             + phizv_*(vtv_*vtv_ + av_*av_*(vxv_*vxv_ + vyv_*vyv_ + vzv_*vzv_)))/av_/av_/(1.0-2.0*phiv_); // replace with stuff
}



/*  ********************************************************  */

/*  This next block of code implements Cash-Karp RKF4-5 based on NRC */

void lightRay::setRKF45values(double h, BaseLens* Lens, int part) {

    if(part == 0){  // used in computing k1
        tv_  = t_;
        xv_  = x_;
        yv_  = y_;
        zv_  = z_;
        vtv_ = vt_;
        vxv_ = vx_;
        vyv_ = vy_;
        vzv_ = vz_;
    } // closes if for RK initial step

    if(part == 1){  // used in computing k2
        tv_  = t_ + b21*k1t_;
        xv_  = x_  + b21*k1x_;
        yv_  = y_  + b21*k1y_;
        zv_  = z_  + b21*k1z_;
        vtv_ = vt_ + b21*k1vt_;
        vxv_ = vx_ + b21*k1vx_;
        vyv_ = vy_ + b21*k1vy_;
        vzv_ = vz_ + b21*k1vz_;
    } // closes if for RK 1st step

    if(part == 2){  // used in computing k3
        tv_  = t_  + b31*k1t_ + b32*k2t_;
        xv_  = x_  + b31*k1x_ + b32*k2x_;
        yv_  = y_  + b31*k1y_ + b32*k2y_;
        zv_  = z_  + b31*k1z_ + b32*k2z_;
        vtv_ = vt_ + b31*k1vt_ + b32*k2vt_;
        vxv_ = vx_ + b31*k1vx_ + b32*k2vx_;
        vyv_ = vy_ + b31*k1vy_ + b32*k2vy_;
        vzv_ = vz_ + b31*k1vz_ + b32*k2vz_;
    } // closes if for RK  2nd step

    if(part == 3){  // used in computing k4
        tv_  = t_  + b41*k1t_ + b42*k2t_ + b43*k3t_;
        xv_  = x_  + b41*k1x_ + b42*k2x_ + b43*k3x_;
        yv_  = y_  + b41*k1y_ + b42*k2y_ + b43*k3y_;
        zv_  = z_  + b41*k1z_ + b42*k2z_ + b43*k3z_;
        vtv_ = vt_ + b41*k1vt_ + b42*k2vt_ + b43*k3vt_;
        vxv_ = vx_ + b41*k1vx_ + b42*k2vx_ + b43*k3vx_;
        vyv_ = vy_ + b41*k1vy_ + b42*k2vy_ + b43*k3vy_;
        vzv_ = vz_ + b41*k1vz_ + b42*k2vz_ + b43*k3vz_;
    } // closes if for RK  3rd step

    if(part == 4){  // used in computing k5
        tv_  = t_  + b51*k1t_ + b52*k2t_ + b53*k3t_ + b54*k4t_;
        xv_  = x_  + b51*k1x_ + b52*k2x_ + b53*k3x_ + b54*k4x_;
        yv_  = y_  + b51*k1y_ + b52*k2y_ + b53*k3y_ + b54*k4y_;
        zv_  = z_  + b51*k1z_ + b52*k2z_ + b53*k3z_ + b54*k4z_;
        vtv_ = vt_ + b51*k1vt_ + b52*k2vt_ + b53*k3vt_ + b54*k4vt_;
        vxv_ = vx_ + b51*k1vx_ + b52*k2vx_ + b53*k3vx_ + b54*k4vx_;
        vyv_ = vy_ + b51*k1vy_ + b52*k2vy_ + b53*k3vy_ + b54*k4vy_;
        vzv_ = vz_ + b51*k1vz_ + b52*k2vz_ + b53*k3vz_ + b54*k4vz_;
    } // closes if for RK  4th step

    if(part == 5){  // used in computing k6
        tv_  = t_  + b61*k1t_ + b62*k2t_ + b63*k3t_ + b64*k4t_ + b65*k5t_;
        xv_  = x_  + b61*k1x_ + b62*k2x_ + b63*k3x_ + b64*k4x_ + b65*k5x_;
        yv_  = y_  + b61*k1y_ + b62*k2y_ + b63*k3y_ + b64*k4y_ + b65*k5y_;
        zv_  = z_  + b61*k1z_ + b62*k2z_ + b63*k3z_ + b64*k4z_ + b65*k5z_;
        vtv_ = vt_ + b61*k1vt_ + b62*k2vt_ + b63*k3vt_ + b64*k4vt_ + b65*k5vt_;
        vxv_ = vx_ + b61*k1vx_ + b62*k2vx_ + b63*k3vx_ + b64*k4vx_ + b65*k5vx_;
        vyv_ = vy_ + b61*k1vy_ + b62*k2vy_ + b63*k3vy_ + b64*k4vy_ + b65*k5vy_;
        vzv_ = vz_ + b61*k1vz_ + b62*k2vz_ + b63*k3vz_ + b64*k4vz_ + b65*k5vz_;
    } // closes if for RK  5th step

    av_ = 1.0;   // these two functions will eventually depend on the tv_ coordinate
    apv_ = 0.0;
    rv_ = sqrt(xv_*xv_ + yv_*yv_ + zv_*zv_);

    phiv_  = Lens->getU(xv_, yv_, zv_);
    phixv_ = Lens->getUdx(xv_, yv_, zv_, h/10.0); // since RKF45, numerical derivative using h/10
    phiyv_ = Lens->getUdy(xv_, yv_, zv_, h/10.0);
    phizv_ = Lens->getUdz(xv_, yv_, zv_, h/10.0);

} // closes function to set values

/* This file takes 1 RKF step */

void lightRay::takeRKF45step(double h, BaseLens* Lens) {

    setRKF45values(h, Lens, 0); // setting up initial values

    k1t_ = h*vtv_;
    k1x_ = h*vxv_;
    k1y_ = h*vyv_;
    k1z_ = h*vzv_;
    k1vt_ = h*computeVtdot();
    k1vx_ = h*computeVxdot();
    k1vy_ = h*computeVydot();
    k1vz_ = h*computeVzdot();

    setRKF45values(h, Lens, 1); // set up values to use in k2's

    k2t_ = h*vtv_;
    k2x_ = h*vxv_;
    k2y_ = h*vyv_;
    k2z_ = h*vzv_;
    k2vt_ = h*computeVtdot();
    k2vx_ = h*computeVxdot();
    k2vy_ = h*computeVydot();
    k2vz_ = h*computeVzdot();

    setRKF45values(h, Lens, 2); // set up values to use in k3's

    k3t_ = h*vtv_;
    k3x_ = h*vxv_;
    k3y_ = h*vyv_;
    k3z_ = h*vzv_;
    k3vt_ = h*computeVtdot();
    k3vx_ = h*computeVxdot();
    k3vy_ = h*computeVydot();
    k3vz_ = h*computeVzdot();

    setRKF45values(h, Lens, 3); // set up values to use in k4's

    k4t_ = h*vtv_;
    k4x_ = h*vxv_;
    k4y_ = h*vyv_;
    k4z_ = h*vzv_;
    k4vt_ = h*computeVtdot();
    k4vx_ = h*computeVxdot();
    k4vy_ = h*computeVydot();
    k4vz_ = h*computeVzdot();

    setRKF45values(h, Lens, 4); // set up values to use in k5's

    k5t_ = h*vtv_;
    k5x_ = h*vxv_;
    k5y_ = h*vyv_;
    k5z_ = h*vzv_;
    k5vt_ = h*computeVtdot();
    k5vx_ = h*computeVxdot();
    k5vy_ = h*computeVydot();
    k5vz_ = h*computeVzdot();

    setRKF45values(h, Lens, 5); // set up values to use in k6's

    k6t_ = h*vtv_;
    k6x_ = h*vxv_;
    k6y_ = h*vyv_;
    k6z_ = h*vzv_;
    k6vt_ = h*computeVtdot();
    k6vx_ = h*computeVxdot();
    k6vy_ = h*computeVydot();
    k6vz_ = h*computeVzdot();

    // save values into temporary spots

    tTemp_ = t_ + c1*k1t_ + c2*k2t_ + c3*k3t_ + c4*k4t_ + c5*k5t_ + c6*k6t_;
    xTemp_ = x_ + c1*k1x_ + c2*k2x_ + c3*k3x_ + c4*k4x_ + c5*k5x_ + c6*k6x_;
    yTemp_ = y_ + c1*k1y_ + c2*k2y_ + c3*k3y_ + c4*k4y_ + c5*k5y_ + c6*k6y_;
    zTemp_ = z_ + c1*k1z_ + c2*k2z_ + c3*k3z_ + c4*k4z_ + c5*k5z_ + c6*k6z_;
    vtTemp_ = vt_ + c1*k1vt_ + c2*k2vt_ + c3*k3vt_ + c4*k4vt_ + c5*k5vt_ + c6*k6vt_;
    vxTemp_ = vx_ + c1*k1vx_ + c2*k2vx_ + c3*k3vx_ + c4*k4vx_ + c5*k5vx_ + c6*k6vx_;
    vyTemp_ = vy_ + c1*k1vy_ + c2*k2vy_ + c3*k3vy_ + c4*k4vy_ + c5*k5vy_ + c6*k6vy_;
    vzTemp_ = vz_ + c1*k1vz_ + c2*k2vz_ + c3*k3vz_ + c4*k4vz_ + c5*k5vz_ + c6*k6vz_;

    // store err values

    errvec[0] = fabs((c1-c1s)*k1t_ + (c2-c2s)*k2t_ + (c3-c3s)*k3t_ + (c4-c4s)*k4t_ + (c5-c5s)*k5t_ + (c6-c6s)*k6t_);
    errvec[1] = fabs((c1-c1s)*k1x_ + (c2-c2s)*k2x_ + (c3-c3s)*k3x_ + (c4-c4s)*k4x_ + (c5-c5s)*k5x_ + (c6-c6s)*k6x_);
    errvec[2] = fabs((c1-c1s)*k1y_ + (c2-c2s)*k2y_ + (c3-c3s)*k3y_ + (c4-c4s)*k4y_ + (c5-c5s)*k5y_ + (c6-c6s)*k6y_);
    errvec[3] = fabs((c1-c1s)*k1z_ + (c2-c2s)*k2z_ + (c3-c3s)*k3z_ + (c4-c4s)*k4z_ + (c5-c5s)*k5z_ + (c6-c6s)*k6z_);
    errvec[4] = fabs((c1-c1s)*k1vt_ + (c2-c2s)*k2vt_ + (c3-c3s)*k3vt_ + (c4-c4s)*k4vt_ + (c5-c5s)*k5vt_ + (c6-c6s)*k6vt_);
    errvec[5] = fabs((c1-c1s)*k1vx_ + (c2-c2s)*k2vx_ + (c3-c3s)*k3vx_ + (c4-c4s)*k4vx_ + (c5-c5s)*k5vx_ + (c6-c6s)*k6vx_);
    errvec[6] = fabs((c1-c1s)*k1vy_ + (c2-c2s)*k2vy_ + (c3-c3s)*k3vy_ + (c4-c4s)*k4vy_ + (c5-c5s)*k5vy_ + (c6-c6s)*k6vy_);
    errvec[7] = fabs((c1-c1s)*k1vz_ + (c2-c2s)*k2vz_ + (c3-c3s)*k3vz_ + (c4-c4s)*k4vz_ + (c5-c5s)*k5vz_ + (c6-c6s)*k6vz_);

}