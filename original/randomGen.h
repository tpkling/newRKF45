
//
// Created by srowell on 4/18/2017.
// Material added by tpk on 4/19/2017 to bring in random galaxy (source) locations
//
#include <random>
#include <time.h>

#ifndef NEW_RKF45_WITH_LENS_RANDOMGEN_H
#define NEW_RKF45_WITH_LENS_RANDOMGEN_H


class randomGen {
public:
    static randomGen& getInstance(){
        static randomGen randomInstance;

        return randomInstance;
    }

    randomGen(randomGen const&) = delete;
    void operator=(randomGen const&) = delete;

    double getRandomRs();
    double getRandomGM();
    double getRandomC();
    double getRandomXpos();
    double getRandomYpos();
    double getRandomZpos();
    double getRandomSourceX();
    double getRandomSourceY();
    double getRandomSourceRedshift();
    void setScale(double scale);
    double getScale();



private:

    double _tScale; // in meters
    double _rsScale;
    double _GMScale; // dimensionless

    std::default_random_engine _re;
    std::normal_distribution<double> _rsLensDistribution;
    std::normal_distribution<double> _GMLensDistribution;
    std::normal_distribution<double> _cLensDistribution;
    std::uniform_real_distribution<double> _xpLensDistribution;
    std::uniform_real_distribution<double> _ypLensDistribution;
    std::uniform_real_distribution<double> _zpLensDistribution;
    std::normal_distribution<double> _galaxyRedshiftDistribution;  //Redshift normal distribution.
    std::uniform_real_distribution<double> _YSourceDistribution;  // value in arcminutes
    std::uniform_real_distribution<double> _XSourceDistribution;  // value in arcminutes


    randomGen(){

        _re = std::default_random_engine(time(0));
        _rsLensDistribution  = std::normal_distribution<double>(7.714194e18, 0.3 * 7.714194e18);
        _GMLensDistribution  = std::normal_distribution<double>(1.482e14, 0.3*1.482e14);
        _cLensDistribution   = std::normal_distribution<double>(8.0, 0.3*8.0);
        _xpLensDistribution  = std::uniform_real_distribution<double>(-7.714194e18/10.0, 7.714194e18/10.0);
        _ypLensDistribution  = std::uniform_real_distribution<double>(-7.714194e18/10.0, 7.714194e18/10.0);
        _zpLensDistribution  = std::uniform_real_distribution<double>(-7.714194e18/10.0, 7.714194e18/10.0);
        _galaxyRedshiftDistribution = std::normal_distribution<double>(0.4,0.3);   //Redshift normal distribution.
        _YSourceDistribution = std::uniform_real_distribution<double>(-5.0,5.0); // values in arc minutes
        _XSourceDistribution = std::uniform_real_distribution<double>(-5.0,5.0); // values in arc minutes

    };


};


#endif //NEW_RKF45_WITH_LENS_RANDOMGEN_H
