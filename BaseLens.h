//
// Created by SRowell on 4/19/2017.
//

#ifndef NEW_RKF45_WITH_LENS_BASELENS_H
#define NEW_RKF45_WITH_LENS_BASELENS_H


class BaseLens {


public:

    /*
     * Abstract class for Lenses
     * All real lenses will override these functions..
     */
    virtual double getU(double, double, double)=0;
    virtual double getUdx(double, double, double, double)=0;
    virtual double getUdy(double, double, double, double)=0;
    virtual double getUdz(double, double, double, double)=0;
    virtual void getPos()=0;

};


#endif //NEW_RKF45_WITH_LENS_BASELENS_H
