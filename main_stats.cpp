#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "modis_hdf.h"
#include "modis_process.h"
#include "alert.h"
#include "testime.h"


extern long converttime_unix(int year, int mon, int day, int hr, int min, int sec) ;
// this program creates a multibanded image centered at the user specified scene center. We then go through each 
// of the MODIS rad and geo files, then based upon b221 value, find the nearest pixel in the grid, add that dn to the total band
// increment that pixels count 

int main(int argc, char *argv[]) {

    char flist[420], infile [420], cfile[420], *tmpname, prefix[420], statsfile[420];
    char flocfile [420];
    char monlog[420], flayers[420];
    char infile_geo [420], outstr[240];
    
    //float *statsdata;
    double *statsdata;
    unsigned short *countdata ;

    float filecount = 0, startlat, startlon, gspace;
    float *lat, *lon, *b21, *b31, *b22, *b32, latval, lonval, xdist, ydist, dist;
    int i, uval, thismonth, ns_modis, nl_modis, npix_modis, datearr[5], modisflag;
    int nx, ny, npix_grid, gridnum;
    long unxtime;
    float cent_lat, cent_lon, alvalue ;
    float b32temp, b22temp ;
    
    alert *al = new alert() ;

    int mon_days[] = {1, 32, 62, 92, 123, 153, 184, 214, 245, 275, 305, 336};
    // count hist is the monthly histogram for each alert
    

    int this_mday;

    if (argc < 3) {
        cout << "Usage: modis_stats flist  log_prefix modisflag c_lon c_lat" << endl;
        exit(-1);
    }


    strcpy(flist, *++argv);
    strcpy(prefix, *++argv);
    modisflag = atoi(*++argv);
    cent_lon = atof(*++argv);
    cent_lat = atof(*++argv);
    
    strcpy(statsfile, prefix);
    strcpy(cfile, prefix);
    strcat(cfile, "_count");

    modis_hdf *therm, *geo;

    float *mnstdv;
    


    // input file is an ascii file, each line has the radiance and the
    // geom file
    // input file contains the full path names of each file 

    nx = 2048;
    ny = 2048;
    ns_modis = 1354;
    nl_modis = 2030;
    // hawaii
    //cent_lat = 19.4 ;
    //cent_lon = -155.25 ;
    // iceland
    //cent_lat = 64.64 ;
    //cent_lon = -17.528 ;
    npix_modis = ns_modis * nl_modis;
    npix_grid = nx * ny ;
    mnstdv = new float [npix_modis * 2];

    // initialize stats data
    countdata = new unsigned short [npix_grid ];
    statsdata = new double [npix_grid* 5];
    
    for (i=0; i<nx * ny * 5; i++) {
        statsdata[i] = 0. ;
        if (i<nx * ny)
            countdata[i]=0 ;
        if (i >= 2*npix_grid && i<3*npix_grid)
            statsdata[i] = 100. ;
        if (i>=4*npix_grid) statsdata[i]=-2. ;
        
    }
    



    FILE *ffile = fopen(flist, "r");


    gspace = .0092;
    startlat = cent_lat + ny / 2. * gspace;
    startlon = cent_lon - nx / 2 * gspace;
    
   

    /*
    sinu_1km *sinu = new sinu_1km();
    sinu->makegrid(basegrid, startlat, startlon, gspace, nx, ny);
    for (i = 0; i < npix_grid; i++) {
        layers[i] = basegrid[i * 2]*100;
        uval = 1 * layers[i];
        if (uval < 0) uval = 0;
        if (uval > 100) uval = 100;
        ucval = (unsigned char) uval;
        layers[npix_grid + i] = basegrid[i * 2 + 1]*1000;
        locdata[i] = ucval;
        locdata[npix_grid + i] = ucval;
        locdata[2 * npix_grid + i] = ucval;
        layers[2 * npix_grid + i] = 0;
        layers[3 * npix_grid + i] = 0;
        layers[4 * npix_grid + i] = 0;
        layers[5 * npix_grid + i] = 0;
    }
     */

    while (!feof(ffile)) {
        fscanf(ffile, "%s %s", infile, infile_geo);
        tmpname = strstr(infile, "021KM") - 3;

        therm = new modis_hdf(infile, modisflag);
        therm->get_date_period(infile, datearr);
        unxtime = converttime_unix(datearr[0], datearr[1], datearr[2], datearr[3], datearr[4], 0);
        al->set_max(therm->b21max, therm->b22max);
        thismonth = datearr[1] - 1;
        this_mday = mon_days[thismonth];
        if (this_mday  > 1) {
            delete therm ;
            break ;
        }
        geo = new modis_hdf(infile_geo, modisflag);
        geo->get_date_period(infile_geo, datearr);
        geo->init_MOD03();

        if (geo->dayflag) {
            cout << "infile  : " << infile << endl << "DAYFILE" << endl;
            delete therm;
            delete geo;
            continue;
        }

        if (!geo->geom_status) {
            //cout << infile_geo << "UHOHOOH"<< endl ;
            cout << "problem with MOD03 FILE" << endl;
        }

        filecount++;
        //if (filecount >10) break ;


        lat = geo->latarr;
        lon = geo->lonarr;
        b21 = &therm->raddata_cal[0];
        b22 = &therm->raddata_cal[npix_modis];
        b32 = &therm->raddata_cal[3 * npix_modis];


        // now fill the layers grids 
        int index, grid_x, grid_y, iline, isamp;
        
        float b221;
        int count = 0, grid_loc;
        
        for (i = 0; i < npix_modis ; i++) {
            // get the image pixel lat and lon
            latval = lat[i];
            lonval = lon[i];
            // find corresponding stats location 
            grid_x = int((lonval - startlon) / gspace + 0.5);
            grid_y = int((startlat - latval) / gspace + 0.5);
            if (grid_x < 0 || grid_y < 0) continue;
            if (grid_x >= nx || grid_y >= ny) continue;
            grid_loc = grid_y * nx + grid_x ;
            
            
            
            b221 = b22[i] ;
            if (b221 >2.) b221 = b21[i] ;
            if (b221 <0.04) continue ;
            alvalue = (b221 - b32[i]) / (b221 + b32[i]) ;
            if (alvalue > -.8) continue;
            b32temp = al->bb_radtotemp (12, b32[i]) ;
            b22temp = al->bb_radtotemp (3.95, b221) ;
            
            
            // update stats array
            countdata[grid_loc]++ ;
            statsdata[grid_loc]+= b221 ;
            statsdata[grid_loc+npix_grid]+=(b221 * b221) ;
            
            
            if (b221 < statsdata[grid_loc+2*npix_grid]){
                statsdata[grid_loc+2*npix_grid] = b221 ;
            }
            if (b22temp > b32temp+30.) continue ;
            if (b221 > statsdata[grid_loc+3*npix_grid]){
                statsdata[grid_loc+3*npix_grid] =b221;
            }
            if (alvalue > statsdata[grid_loc+4*npix_grid]){
                statsdata[grid_loc+4*npix_grid] =alvalue;
            }
            

        }
        
        
        
        delete therm;
        delete geo;
    } // read in the next pair of modis scenes

    
      // now calc the mean h
    int N;
    double meanval, stdval ;
      for (i=0;i<npix_grid; i++) {
            if (countdata[i] > 0.) {
                
                N = countdata[i] ;
                meanval = statsdata[i] / N ;
                stdval = sqrt( (statsdata[i+npix_grid] - (statsdata[i] * statsdata[i]/N) )/N ) ;
                statsdata[i] = meanval ;
                statsdata[i+npix_grid]=stdval ;
            }
      }
    
    
    FILE *fout = fopen(statsfile, "w") ;
    fwrite (statsdata, 8, npix_grid * 5, fout) ;
    fclose (fout) ;
    FILE *foutc = fopen(cfile, "w") ;
    fwrite (countdata, 2, npix_grid , fout) ;
    fclose (fout) ;

	
    
    char ohdr[240];
    strcpy(ohdr, statsfile);
    strcat(ohdr, ".hdr");
    al->write_envi_header(ohdr, nx, ny, startlat, startlon, gspace);
    


    delete [] statsdata ;
    delete [] countdata ;
    delete al ;

}
