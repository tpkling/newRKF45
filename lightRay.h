//
// Developed by Joey Garuti, Tyler Martell, Amanda Coughlin, Anne Marie Pooler & Adam St. Amand
// Based on tpk code for PHYS 422
// Modified to work with lens constructor from PHYS 422 by tpk on April 19, 2017
//

#ifndef NEW_RKF45_WITH_LENS_LIGHTRAY_H
#define NEW_RKF45_WITH_LENS_LIGHTRAY_H

#include "BaseLens.h"
#include "Lens.h"

class lightRay {

private:

public:
    // define array for error

    double errvec[8];
    double errmax = 0.0;


    // define the coordinates

    double t_, tTemp_;
    double x_, xTemp_;
    double y_, yTemp_;
    double z_, zTemp_;
    double vt_, vtTemp_;
    double vx_, vxTemp_;
    double vy_, vyTemp_;
    double vz_, vzTemp_;

    // define RK4 k's

    double k1t_, k2t_, k3t_, k4t_, k5t_, k6t_;
    double k1x_, k2x_, k3x_, k4x_, k5x_, k6x_;
    double k1y_, k2y_, k3y_, k4y_, k5y_, k6y_;
    double k1z_, k2z_, k3z_, k4z_, k5z_, k6z_;
    double k1vt_, k2vt_, k3vt_, k4vt_, k5vt_, k6vt_;
    double k1vx_, k2vx_, k3vx_, k4vx_, k5vx_, k6vx_;
    double k1vy_, k2vy_, k3vy_, k4vy_, k5vy_, k6vy_;
    double k1vz_, k2vz_, k3vz_, k4vz_, k5vz_, k6vz_;

    // define variables for calculations

    double tv_, xv_, yv_, zv_, vtv_, vxv_, vyv_, vzv_;
    double phiv_, phixv_, phiyv_, phizv_;
    double rv_, av_, apv_; // for cosmological expansion

    // define the constructors

    lightRay(); // default constructor
    lightRay(double t, double x, double y, double z, double vt, double vx, double vy, double vz);

    // define functions to compute terms in ODEs
    // these use the public variable data, so we don't need to pass in any variables

    double computeVtdot();
    double computeVxdot();
    double computeVydot();
    double computeVzdot();

    // a function to set the middle values in the RKF 4-5 code

//    void setRKF45values(double h, double sigmav_, double rc, int part);
    void setRKF45values (double h, BaseLens* oneLens, int part);
    // a function to take a Runge Kutta Fehlberg step

//    void takeRKF45step(double h, double sigmaV, double rc);
    void takeRKF45step(double h, BaseLens* oneLens);

};

#endif //NEW_RKF45_WITH_LENS_LIGHTRAY_H
