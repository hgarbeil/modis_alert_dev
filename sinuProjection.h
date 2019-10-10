/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinuProjection.h
 * Purpose : Header file for the sinusoidal projection class. This class is used to relate latitude longitude to an xy 
 * system based upon the sinusoidal projection, commonly used for modis global daasets.
 * Author: harold garbeil
 *
 * Created on August 3, 2018, 9:29 AM
 */

#ifndef SINUPROJECTION_H
#define SINUPROJECTION_H

#include <math.h>
#include <iostream>
#include <vector>

using namespace std ;

class sinuProjection {
public:
    sinuProjection();
    sinuProjection(const sinuProjection& orig);
    virtual ~sinuProjection();
    void set_parameters (float radius, float centmeridian)  ;
    void fillGrid (int mon, int aq_terra_flag, float ulc_lat, float ulc_lon, float gridspace, int nx, int ny, float *outarr) ; 
    void latlon_to_xy (float lat, float lon, float *x, float *y) ;
    void xy_to_latlon (float x, float y, float *lat, float *lon) ;
    //void getValues(float lat, float lon, float *vals) ;
    void getCorrespondingValues(vector<int> alinds, int mon, int aq_terra_flag, float *lat, float *lon, float *vals) ;
    void getNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *outarr, bool openFlag) ;
    void getAllValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outNTI, bool openFlag) ;
    void getMaxNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outmaxNTI, bool openFlag) ;
    float mDegToRad, wrld_startx, wrld_starty, gspace ;
    char basef[420] ;
private:
    float Radius ;
    float cent_merid ;
};

#endif /* SINUPROJECTION_H */

