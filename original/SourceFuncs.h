//
// Created by TKLING on 4/19/2017.
//

#ifndef NEW_RKF45_WITH_LENS_SOURCEFUNCS_H
#define NEW_RKF45_WITH_LENS_SOURCEFUNCS_H


double computeLookbackTime(double zr);
double a(double t);
double zdot(double t);
double zPositionIntegrator(double t1, double t2, double z1);  // to integrate from t1 to t2 using standard RK 4th order


#endif //NEW_RKF45_WITH_LENS_SOURCEFUNCS_H
