#include <iostream>
#include "SourceFuncs.h"
#include "Lens.h"
#include "LensGroup.h"
#include "BaseLens.h"
#include "randomGen.h"
#include "lightRay.h"
#include <cmath>

using namespace std;

int main() {

    // Code snippets following this are modified from work of Maria and Nick

    //create Singleton used for random number generating..
    // randomGen::getInstance();

    double xSource = randomGen::getInstance().getRandomSourceX();
    double ySource = randomGen::getInstance().getRandomSourceY();  // in arc minutes
    double redshiftLens = 0.45;

    double redshiftSource = randomGen::getInstance().getRandomSourceRedshift();

    while(redshiftSource < 0.6 || redshiftSource> 1.0) {
        redshiftSource = randomGen::getInstance().getRandomSourceRedshift();
    }

    double tNow = computeLookbackTime(0.0); // value of now time in meters
    double tLens = computeLookbackTime(redshiftLens); // value of lens time in meters
    double tSource = computeLookbackTime(redshiftSource); // value of source time in meters

    double zobs = -1.0*zPositionIntegrator(tNow, tLens, 0.0); // places observer on + z axis
    double zSource = zPositionIntegrator(tLens, tSource, 0.0); // places source in -z position

    double tScale = tNow;
    double zObsScale = zobs/tScale;
    double zSourceScale = zSource/tScale;

    double xSourceScale = xSource*(-zSourceScale+zObsScale)*0.00029088; // number converts arc minutes to radians
    double ySourceScale = ySource*(-zSourceScale+zObsScale)*0.00029088; // number converts arc minutes to radians

    cout<< xSourceScale<<"  "<< ySourceScale <<"  "<< zSourceScale <<"  "<< zObsScale<<endl;

    randomGen::getInstance().setScale(tScale);
    cout<<randomGen::getInstance().getScale()<<endl;


    //
    // This test code was written by Steven Rowell
    // Example code illustrating the use of a lens cluster and a single lens and storing them inside a container of the parent class.
    //

    //Two pointers to a BaseLens object
    BaseLens *clusterLens;
    BaseLens *oneLens;



    //Implementing a single lens and lens group!
    Lens singleLens = Lens();
    LensGroup lotsOfLens = LensGroup(5);  // the number here in the () indicates number of clumps

    //Assigning references to the pointers
    oneLens = &singleLens;
    clusterLens = &lotsOfLens;

    //Creating a container for our lenses
    //BaseLens * set[] = {&lotsOfLens,&singleLens}; would also work..
    BaseLens * set[] = {clusterLens, oneLens};

    std::cout <<  clusterLens->getUdx(0.2,0.1,0.2, 1.0e-3) << std::endl;
    std::cout << (*set[0]).getUdx(0.2,0.1,0.2, 1.0e-3)  << std::endl;

    std::cout <<  oneLens->getUdx(0.2,0.1,0.2, 1.0e-3)  << std::endl;
    std::cout << (*set[1]).getUdx(0.2,0.1,0.2, 1.0e-3)  << std::endl;

    //Uncomment to dump lens data to console..
    (*set[0]).getPos();

// The code below here was worked on by the RKF group

    double t = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = zObsScale;
    double vt = 1.0;
    double vy = 0.0;
    double vx = sin(0.03);
    double vz = -cos(0.03); // since initial z > 0, want to come towards the origin with negative z values

    double errmax = 0.0;
    double err = 1e-2;

    lightRay ray(t, x, y, z, vt, vx, vy, vz);

    double hnew;
    double h = 0.001;
    double rc = 1.0;
    double sigmaV = 0.5;

    int countSafe = 0;

    while (ray.z_ > -20.0) {
        countSafe++;
        errmax = 0.0;

        ray.takeRKF45step(h, clusterLens); // runs a ray for one step

        //std::cout << "intial h = " << h << std::endl;

        for (int i = 0; i < 8; i++) {
            if (ray.errvec[i] > errmax) {
                errmax = ray.errvec[i];
                //errmax = 1e-4;
                //std::cout << "errmax = " << errmax << std::endl;
            }
        }

        if (err / errmax > 1.0) {
            ray.t_ = ray.tTemp_;
            ray.x_ = ray.xTemp_;
            ray.y_ = ray.yTemp_;
            ray.z_ = ray.zTemp_;
            ray.vt_ = ray.vtTemp_;
            ray.vx_ = ray.vxTemp_;
            ray.vy_ = ray.vyTemp_;
            ray.vz_ = ray.vzTemp_;
            hnew = h*pow((err/errmax), 0.2);
            if (hnew>1.3*h) {
                h = 1.3*h;
            }
            else h = hnew;

            //std::cout << "h get bigger = " << h << std::endl;
        }
        else {
            hnew = h*pow(err/errmax, 0.2);
            if(hnew < 0.7*h){
                h = 0.7*h;
            }
            else h = hnew;
            //std::cout << "h get smaller = " << h << std::endl;
        }

        if (countSafe % 4 == 0) {
            // std::cout << ray.errvec[0] << " " << ray.errvec[1] << " " << ray.errvec[2] << " " << ray.errvec[3] << " "
            //<< ray.errvec[4] << " " << ray.errvec[5] << " " << ray.errvec[6] << " " << ray.errvec[7] << std::endl;
            std::cout << ray.t_ << "  " << ray.x_ << "  " << ray.y_ << "  " << ray.z_ << std::endl;
        }
        if(countSafe > 20) break;
    }//closes while loop

    return 0;
}