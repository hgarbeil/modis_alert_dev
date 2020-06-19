#ifndef mhdf
#define mhdf 

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include "mfhdf.h"

using namespace std;

class modis_hdf {
    char inhdfname [420];
    float *modbands_em, *modbands_refl;

    int32 sd_id, n_datasets, n_fileattrs;


public:
    bool thermFlag, aquaflag, dayflag, geom_status ;
    int aq_terra_flag, readstatus, modflag, hour ;
	unsigned char *badpix ;
    uint16 *thermdata, *refldata;
    float *latarr, *lonarr;
    int16 *solzen, *solaz, *senszen, *sensaz;
    float *th_scales_offsets, *refl_scales_offsets;
    float *raddata_cal, *refdata_cal, *solsens;
	float *glint ;
	float b21max, b22max ;
    modis_hdf(char *, int);
    ~modis_hdf();
    int open_hdf();
    void get_aq_terra_flag () ;
    void init_MOD21();
    void init_MOD03();
	int  calc_glint() ;
    int load_geometry();
    int load_thermal_bands();
    void calib_thermal_bands();
    void load_refSB_bands();
    void calib_refSB_bands();
    void calib_geometry();
    void get_date_period(char *infile, int *datearr);
    void get_date_period_fix(char *infile, int *datearr);
    void get_file_name(char *, char *);



};




#endif

	

