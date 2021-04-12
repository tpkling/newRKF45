//
// File created by TKLING on 4/19/2017.
// Code written by Maria Patrone and Nick Matsuo in PHYS 422
//
#include<cmath>
#include<iostream>

double computeLookbackTime(double zr){
    double omegav = 0.7;
    double omegam = 0.3;
    double a = 1.0/(1.0 + zr);
    double Ho = 2.25e-18;

    return (2.0 * asinh(pow(pow(omegav/omegam, 1.0/3.0) * a, 1.5)))/(3.0 * Ho * sqrt(omegav));
}

double a(double t){
    double omegam = 0.3;
    double omegav = 0.7;
    double Ho = 2.25e-18;

    return pow(omegam/omegav, 1.0/3.0) * pow(sinh((3.0 * Ho * sqrt(omegav) * t)/2.0), 2.0/3.0);
}

double zdot(double t){
    return 1.0/a(t);
}

double zPositionIntegrator (double t1, double t2, double z1){

    double h = -0.01 * (t1 - t2);  // t2 < t1 ::  t2 is more in the past
    double k1, k2, k3, k4, t, z;

    z = z1;
    t = t1;

    for(int i=0; i<100; i++){

        k1 = zdot(t);
        k2 = zdot(t + 0.5 * h);
        k3 = zdot(t + 0.5 * h);
        k4 = zdot(t + h);

        z = z + h/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t = t + h;
        //std::cout<< t << "  " << z << std::endl;

    }

    std::cout<< "Did we get to the desired time? " << t2 <<"  "<< t<<"  "<< (t2-t)/t2<< std::endl;
    return z;
}