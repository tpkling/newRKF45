//
// Created by srowell on 4/18/2017.
//

#include "randomGen.h"
#include<iostream>

double randomGen::getRandomRs() {
    return _rsLensDistribution(_re);

}

double randomGen::getRandomGM() {
    return _GMLensDistribution(_re);
}

double randomGen::getRandomC() {
    return _cLensDistribution(_re);
}

double randomGen::getRandomXpos() {
    return _xpLensDistribution(_re);
}

double randomGen::getRandomYpos() {
    return _ypLensDistribution(_re);
}

double randomGen::getRandomZpos() {
    return _zpLensDistribution(_re);
}

double randomGen::getRandomSourceX()  {
    return _XSourceDistribution(_re);
}

double randomGen::getRandomSourceY()  {
    return _YSourceDistribution(_re);
}

double randomGen::getRandomSourceRedshift()  {
    return _galaxyRedshiftDistribution(_re);
}

void randomGen::setScale(double scale) {
    _tScale = scale;
    _GMScale = 1.482e14/scale;
    _rsScale = 7.714194e18/scale;
    std::cout<<_tScale<<std::endl;
}

double randomGen::getScale()  {
    return _tScale;
    std::cout<<_tScale<<std::endl;
}

