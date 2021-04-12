//
// Created by TKLING on 4/19/2017 from the work of students in PHYS 422
// Steven Rowell, Henry Tran, John L'Heureux, Will Coon & Krista Roderick
//

#ifndef NEW_RKF45_WITH_LENS_LENS_H
#define NEW_RKF45_WITH_LENS_LENS_H

#include "BaseLens.h"

class Lens: public BaseLens{
public:
    Lens(); //Default constructor
    Lens(double xp,double yp,double zp, double GM, double rs, double c); //Constructor for lens with known conditions
    double getU(double, double, double);  //Getter for the potential at a point
    double getUdx(double, double, double, double); //Getter for X derivative of potential at a point with a step
    double getUdy(double, double, double, double); //Getter for Y derivative of potential at a point with a step
    double getUdz(double, double, double, double); //Getter for Z derivative of potential at a point with a step
    void getPos();
private:
    double _xp;  //X center
    double _yp;  //Y center
    double _zp;  //Z center
    double _GM;  //const
    double _rs;  //const
    double _c;  //const

};

#endif //NEW_RKF45_WITH_LENS_LENS_H
