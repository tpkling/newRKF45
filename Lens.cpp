//
// Created by TKLING on 4/19/2017 from the work of students in PHYS 422
// Steven Rowell, Henry Tran, John L'Heureux, Will Coon & Krista Roderick
//
//


#include "Lens.h"
#include <cmath>
#include "randomGen.h"
#include <iostream>

static double _PI = acos(-1.0);

Lens::Lens() {
    _xp = randomGen::getInstance().getRandomXpos();
    _yp = randomGen::getInstance().getRandomYpos();
    _zp = randomGen::getInstance().getRandomZpos();
    _GM = randomGen::getInstance().getRandomGM();
    _rs = randomGen::getInstance().getRandomRs();
    _c  = randomGen::getInstance().getRandomC();
}

//TODO Modify constructor for random gen... Also we need the other params added!!!!! Someone will have to deal with this ASAP..
Lens::Lens(double xp,double yp,double zp, double GM, double rs, double c){
    _xp = xp;
    _yp = yp;
    _zp = zp;
    _GM = GM;
    _rs = rs;
    _c = c;

};

double Lens::getU(double x, double y, double z){

    double r = sqrt((x-_xp)*(x-_xp) + (y-_yp)*(y-_yp) + (z-_zp)*(z-_zp));
    //TODO: Not sure why t = 3*_c..
    double t = 3*_c;
    double xR = r/_rs;

    return  _GM/_rs*t*t/(1+t*t)/(1+t*t) *
            (atan(xR/t)*((1/t) -t - (2*t/xR))
             + log((1+xR*xR/t/t)/(1+xR)/(1+xR)) * ((t*t -1)/2/xR - 1)
             + _PI *(t*t -1)/2/t
             - 2*log(t));
}
double Lens::getUdx(double x, double y, double z, double h){
    return (getU(x+h, y, z) - getU(x-h, y, z))/2/h;
}

double Lens::getUdy(double x, double y, double z, double h){
    return (getU(x, y+h, z) - getU(x, y-h, z))/2/h;
}

double Lens::getUdz(double x, double y, double z, double h){
    return (getU(x, y, z+h) - getU(x, y, z-h))/2/h;
}

void Lens::getPos() {
    std::cout << "X " << _xp << " Y " << _yp << " z " << _zp << " GM " << _GM << " RS " << _rs << " C " << _c << std::endl;
}


