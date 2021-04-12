#include <iostream>
#include "BaseLens.h"
#include "Lens.h"
#include "Coords.h"
#include "lightRay.h"
#include <cmath>

using namespace std;

int main() {

    // Modification of the files (one direction down) where we wanted to have more complicated lense
    // This program works on simple truncated NFW at the origin a distance away, no randomization, no lens group

    //Implementing a single lens at (0,0,0) with (GM, rs, c) ! Putting all numbers in default constructor

    Lens singleLens = Lens();
    //Assigning references to the pointers
    BaseLens *oneLens;
    oneLens = &singleLens;

    double c_light = 299792458; // in  m/s
    double t_s = 4.24989e+17; // in s
    t_s = t_s*c_light; // in m
    double kpc    = 3.086e+19; // 1 kpc in meters
    
    cout << 100*kpc/t_s<<endl;
    
    lightRay ray;
    ray.runRay(oneLens);

    return 0;
}

