#ifndef m_rs 
#define m_rs

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "modis_hdf.h"
#include "sinuProjection.h"

using namespace std;

class modis_process {
    bool alertFlag ;
    char modid [44], proclog[420] ;
    
	
	float startlat, startlon, endlat, endlon, gspace ;
	float *distarr, *bandsarr, *temp, *alerts,  *newalerts, *glob_mnstdev, *basedata, *basevals_full, *basevals_full_new ;
	float day_limit, night_limit ;
	modis_hdf *geom, *therm ;
     
	// alert indices
	vector <int> alinds ;


public:
    unsigned short *histoarr ;
    int nsamps, nlines, ns, nl, nx_grid, ny_grid ;
    bool dayFlag ;
    float gridspace;
    modis_process();
    ~modis_process();
    sinuProjection *sp ;
    void write_output(char *);
    void set_month(int m);
    void set_prefix (char *pref) ;
    void set_proclog (char *plog) ;
    void set_modis_hdfs(modis_hdf *geom, modis_hdf *therm);
    void set_bounds(float ulc_lat, float ulc_lon, float llc_lat,
            float llc_lon, float gridspace);
    void set_bounds (float ulc_lat, float ulc_lon,float gridspace, int nl, int ns) ;
    float bb_radtotemp(float, float);
    float bb_temptorad(float, float);
    void process();
    //void procalert();
    void write_header(char *, int);
    void write_alert_textfile(char *);
    void extract_from_baseline_file();
    float calcdist(float, float);
    void calc_alert(char *, char *, bool);
    void calc_alert_test (char *, char *);
    void calc_alert_bb (char *, char *) ;
    void zscore () ;
    void get_nearest_pixel(float lat, float lon, int *ind, float *pargs);
    int mon, num_alerts_full ;
                        


} ;

#endif
