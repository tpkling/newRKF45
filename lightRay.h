//
// Developed by Joey Garuti, Tyler Martell, Amanda Coughlin, Anne Marie Pooler & Adam St. Amand
// Based on tpk code for PHYS 422
// Modified to work with lens constructor from PHYS 422 by tpk on April 19, 2017
//

#ifndef NEW_RKF45_WITH_LENS_LIGHTRAY_H
#define NEW_RKF45_WITH_LENS_LIGHTRAY_H

#include "BaseLens.h"
#include "Lens.h"
#include "Coords.h"

class lightRay {

private:

public:


    Coords coords;

    // define the coordinates

    // define RKF 4-5 k's

    double k1t_, k2t_, k3t_, k4t_, k5t_, k6t_;
    double k1x_, k2x_, k3x_, k4x_, k5x_, k6x_;
    double k1y_, k2y_, k3y_, k4y_, k5y_, k6y_;
    double k1z_, k2z_, k3z_, k4z_, k5z_, k6z_;
    double k1vt_, k2vt_, k3vt_, k4vt_, k5vt_, k6vt_;
    double k1vx_, k2vx_, k3vx_, k4vx_, k5vx_, k6vx_;
    double k1vy_, k2vy_, k3vy_, k4vy_, k5vy_, k6vy_;
    double k1vz_, k2vz_, k3vz_, k4vz_, k5vz_, k6vz_;
    
    double k1Xt_, k2Xt_, k3Xt_, k4Xt_, k5Xt_, k6Xt_;
    double k1Xx_, k2Xx_, k3Xx_, k4Xx_, k5Xx_, k6Xx_;
    double k1Xy_, k2Xy_, k3Xy_, k4Xy_, k5Xy_, k6Xy_;
    double k1Xz_, k2Xz_, k3Xz_, k4Xz_, k5Xz_, k6Xz_;
    double k1VXt_, k2VXt_, k3VXt_, k4VXt_, k5VXt_, k6VXt_;
    double k1VXx_, k2VXx_, k3VXx_, k4VXx_, k5VXx_, k6VXx_;
    double k1VXy_, k2VXy_, k3VXy_, k4VXy_, k5VXy_, k6VXy_;
    double k1VXz_, k2VXz_, k3VXz_, k4VXz_, k5VXz_, k6VXz_;

    // define variables for calculations

    double tv_, xv_, yv_, zv_, vtv_, vxv_, vyv_, vzv_;
    double Xtv_, Xxv_, Xyv_, Xzv_, VXtv_, VXxv_, VXyv_, VXzv_;
    double phiv_, phixv_, phiyv_, phizv_; // value of phi and the 1st derivatives
    double phixxv_, phiyyv_, phizzv_, phixyv_, phixzv_, phiyzv_; // all the second derivatives
    double R0301_, R0302_, R0303_, R0331_, R0332_; // t Riemann tenson terms needed in code
    double R1001_, R1331_, R1031_, R1002_, R1030_, R1301_, R1330_, R1332_; // x
    double R2002_, R2332_, R2032_, R2001_, R2030_, R2302_, R2330_, R2331_; // y
    double R3030_, R3032_, R3031_, R3002_, R3001_; // z
    
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
    
    double computeVXtdot();
    double computeVXxdot();
    double computeVXydot();
    double computeVXzdot();

    // a function to set the middle values in the RKF 4-5 code
    void setRKF45values (Coords coords, double h, BaseLens* oneLens, int part);
    // a function to take a Runge Kutta Fehlberg step
    Coords takeRKF45step(Coords coords, double errvec[], double h, BaseLens* oneLens);  // returns the new coords   
    void runRay(BaseLens* oneLens);

};

#endif //NEW_RKF45_WITH_LENS_LIGHTRAY_H
