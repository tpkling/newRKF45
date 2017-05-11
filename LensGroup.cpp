//
// Created by srowell on 4/11/2017.
//

#include "LensGroup.h"
#include <iostream>


LensGroup::LensGroup(int n) {

    //TODO: Randomly initialize parameters of the created lenses

    _GroupSize = n;          //store how many lenses are contained in group

    _LensSet = new Lens[n];  //Init a vector of n lenses

    for (int m = 0; m < n; m++) {
        _LensSet[m] = Lens();
    }

}

double LensGroup::getU(double x, double y, double z) {

    double lensSum = 0.00;
    for (int m = 0; m < _GroupSize; m++) {
        lensSum+= _LensSet[m].getU(x,y,z);
    }

    return lensSum;
}

double LensGroup::getUdx(double x, double y, double z, double h) {

    double lensSum = 0.00;

    for (int m = 0; m < _GroupSize; m++) {

        lensSum+= _LensSet[m].getUdx(x,y,z,h);
    }

    return lensSum;

}

double LensGroup::getUdy(double x, double y, double z, double h) {

    double lensSum = 0.00;

    for (int m = 0; m < _GroupSize; m++) {
        lensSum+= _LensSet[m].getUdy(x,y,z,h);
    }

    return lensSum;

}

double LensGroup::getUdz(double x, double y, double z, double h) {

    double lensSum = 0.00;

    for (int m = 0; m < _GroupSize; m++) {
        lensSum+= _LensSet[m].getUdz(x,y,z,h);
    }

    return lensSum;

}

void LensGroup::getPos() {
    for (int m = 0; m < _GroupSize; m++) {
        _LensSet[m].getPos();
    }

}





