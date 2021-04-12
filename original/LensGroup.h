//
// Created by srowell on 4/11/2017.
//

#ifndef NEW_RKF45_WITH_LENS_LENSGROUP_H
#define NEW_RKF45_WITH_LENS_LENSGROUP_H

#include "BaseLens.h"
#include "Lens.h"

class LensGroup :  public BaseLens {

public:
    LensGroup(int n);  //constructor creates a set of n lenses

    // these need x, y, z positions and a derivative width - to be set as RKF stepsize / 10 in context of rays
    double getU(double, double, double);            //Getter for potential of n lenses at a point
    double getUdx(double, double, double, double);  //Getter X derivative of set of lenses at a point given a step
    double getUdy(double, double, double, double);  //Getter Y derivative of set of lenses at a point given a step
    double getUdz(double, double, double, double);  //Getter Z derivative of set of lenses at a point given a step

    void getPos(); //Debug function..

private:
    int _GroupSize; //container for lenses
    Lens * _LensSet;



};

#endif //NEW_RKF45_WITH_LENS_LENSGROUP_H
