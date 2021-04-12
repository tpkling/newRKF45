//
// Created by TKLING on 4/19/2017 from the work of students in PHYS 422
// Steven Rowell, Henry Tran, John L'Heureux, Will Coon & Krista Roderick
//
//


#include "Lens.h"
#include <cmath>
#include <iostream>

//static double c1a 0.7539474411 // For H=70, Ol = 0.7
//static double c2a 1.209935122
static double _H   = 70.0/(3.0857e+19)*4.24989e+17; // H/c_light today -- scaled in appropriate units for numerics
static double _PI = acos(-1.0);

Lens::Lens() {
// not using this at the moment, no random generation
//    _xp = randomGen::getInstance().getRandomXpos();
//    _yp = randomGen::getInstance().getRandomYpos();
//    _zp = randomGen::getInstance().getRandomZpos();
//    _GM = randomGen::getInstance().getRandomGM();
//    _rs = randomGen::getInstance().getRandomRs();
//    _c  = randomGen::getInstance().getRandomC();
    
    double c_light = 299792458; // in  m/s
    double t_s = 4.24989e+17; // in s
    t_s = t_s*c_light; // in m
    
//    double Mpc = 3.086e+22; // 1 Mpc in meters
    double kpc    = 3.086e+19; // 1 kpc in meters
    double OmegaM = 0.3;
    double OmegaL = 0.7;
    _av           = 1.0/(1+0.45); // lens at redshift 0.45
    double Hsq    = _H*_H*(OmegaM/pow(_av, 3) + OmegaL);
    _c            = 7.315;
    _rs           = 250*kpc/t_s; // scaled
    double deltaC = 200.0/3.0 * pow(_c, 3)/(log(1+_c) - _c/(1+ _c));
    _GM           = 1.5*deltaC*Hsq*pow(_rs, 3.0); // scaled
    _xp           = 0.0;
    _yp           = 0.0;
    _zp           = 0.0;
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
    double t = 3*_c;
    double xR = _av*r/_rs;

    return  _GM/_rs*t*t/(1+t*t)/(1+t*t) *
            (atan(xR/t)*((1/t) -t - (2*t/xR))
             + log((1+xR*xR/t/t)/(1+xR)/(1+xR)) * ((t*t -1)/2/xR - 1)
             + _PI *(t*t -1)/2/t
             - 2*log(t));
}
double Lens::getUdx(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x+d, y, z) - getU(x-d, y, z))/2/d;
}

double Lens::getUdy(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x, y+d, z) - getU(x, y-d, z))/2/d;
}

double Lens::getUdz(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x, y, z+d) - getU(x, y, z-d))/2/d;
}

double Lens::getUdxx(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x+d, y, z) + getU(x-d, y, z) -2.0*getU(x,y,z) )/pow(d,2.0); 
}

double Lens::getUdyy(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x, y+d, z) + getU(x, y-d, z) -2.0*getU(x,y,z) )/pow(d,2.0); 
}

double Lens::getUdzz(double x, double y, double z, double h){
    double d = h/10.0;
    return (getU(x, y, z+d) + getU(x, y, z-d) -2.0*getU(x,y,z) )/pow(d,2.0); 
}

double Lens::getUdxy(double x, double y, double z, double h){
    double d = h/10.0;
    return 0.25 * (getU(x+d, y+d, z) + getU(x-d, y-d, z) - getU(x+d,y-d,z) - getU(x-d,y+d,z) )/pow(d,2.0); 
}

double Lens::getUdxz(double x, double y, double z, double h){
    double d = h/10.0;
    return 0.25 * (getU(x+d, y, z+d) + getU(x-d, y, z-d) - getU(x+d,y,z-d) - getU(x-d,y,z+d) )/pow(d,2.0); 
}

double Lens::getUdyz(double x, double y, double z, double h){
    double d = h/10.0;
    return 0.25 * (getU(x, y+d, z+d) + getU(x, y-d, z-d) - getU(x,y-d,z+d) - getU(x,y+d,z-d) )/pow(d,2.0); 
}


void Lens::getPos() {
    std::cout << "X " << _xp << " Y " << _yp << " z " << _zp << " GM " << _GM << " RS " << _rs << " C " << _c << std::endl;
}


