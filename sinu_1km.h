/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinu_1km.h
 * Purpose : Header file for the sinusoidal projection class. This class is used to relate latitude longitude to an xy 
 * system based upon the sinusoidal projection, commonly used for modis global daasets.
 * Author: harold garbeil
 *
 * Created on August 3, 2018, 9:29 AM
 */

#ifndef SINU1KM_H
#define SINU1KM_H

#include <math.h>
#include <iostream>
#include <vector>

using namespace std ;

class sinu_1km {
public:
    sinu_1km();
    sinu_1km(const sinu_1km& orig);
    virtual ~sinu_1km();
    void set_parameters (float radius, float centmeridian)  ;
	void set_mday (int m) ;
    void fillGrid (int mon, int aq_terra_flag, float ulc_lat, float ulc_lon, float gridspace, int nx, int ny, float *outarr) ; 
    void latlon_to_xy (float lat, float lon, float *x, float *y) ;
    void xy_to_latlon (float x, float y, float *lat, float *lon) ;
    //void getValues(float lat, float lon, float *vals) ;
    void getCorrespondingValues(vector<int> alinds, int mon, int aq_terra_flag, float *lat, float *lon, float *vals) ;
    void getNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *outarr, bool openFlag) ;
    void getAllValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outNTI, bool openFlag) ;
    void getMaxNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outmaxNTI, bool openFlag) ;
    void makegrid (float *basegrid, float startlat, float startlon, float gspace, int nx, int ny) ;
void makegrid_32 (float *basegrid, float startlat, float startlon, float gspace, int nx, int ny) ;
int get_mn_stdev_array (float *lat, float *lon, int inpix, int modisflag, float *mnstddev) ;
int get_mn_stdev_array_old (float *lat, float *lon, int inpix, int modisflag, float *mnstddev) ;
int get_mn_stdev_array_b32 (float *lat, float *lon, int inpix, int modisflag, float *mnstddev) ;
int get_mn_stdev_array_nti (float *lat, float *lon, int inpix, int modisflag, float *mnstddev) ;
int get_gridnum (float lat, float lon, int *x, int *y) ;	
    float mDegToRad, wrld_startx, wrld_starty, gspace ;
	void set_badpix (unsigned char *) ;
	int mday ;
    char basef[420] ;
	unsigned char *badpix ;
private:
    float Radius ;
    float cent_merid, meters_per_degree ;
	float gspace10 ;
};

#endif /* SINUPROJECTION_H */

