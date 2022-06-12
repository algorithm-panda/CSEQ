/* 
 * File:   utility.h
 * Author: hjf
 *
 * Created on April 25, 2016, 10:28 PM
 */

#ifndef UTILITY_H
#define	UTILITY_H

#include<cmath>
#ifdef	__cplusplus
extern "C" {
#endif

#define pi 3.1415926535897932384626433832795

#define EARTH_RADIUS 6378.137 //地球半径 KM

    double rad(double d) {
        return d * pi / 180.0;
    }

    double Spherical_distance(double lat1, double lng1, double lat2, double lng2) {
     //   return sqrt(pow((lat1- lat2),2) + pow((lng1 - lng2),2));

        double radLat1 = rad(lat1);
        double radLat2 = rad(lat2);
        double a = radLat1 - radLat2;
        double b = rad(lng1) - rad(lng2);

        double s = 2 * asin(sqrt(pow(sin(a / 2), 2) + cos(radLat1) * cos(radLat2) * pow(sin(b / 2), 2)));
        s = s * EARTH_RADIUS * 1000;
        return s;
    }


#ifdef	__cplusplus
}
#endif

#endif	/* UTILITY_H */

