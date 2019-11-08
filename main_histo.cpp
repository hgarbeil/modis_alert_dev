#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "modis_hdf.h"
#include "modis_process.h"
#include "alert.h"
#include "testime.h"
#include <vector>

extern int mon_days[] ;

int main (int argc, char *argv[]) {

	char flist[420], infile [420], *tmpname,  prefix[420] ;
	char monlog[420], fhisto[420], fcount[420],fbase[420] ;
	char infile_geo [420], outstr[240]  ;
	
	float filecount =0, startlat, startlon, gspace ;
	float *lat, *lon, *b21, *b32, *b22, latval, lonval, xdist, ydist, dist ;
	float b221val, b32val, Tb221val, Tb32val ;
	int i, countnum, thismonth, ns_modis, nl_modis, npix_modis, datearr[5], modisflag ;
	int lineloc, samploc, newloc ;
	int npix_grid, gridnum ;
	long unxtime ;
	float ulc_lat, ulc_lon, lrc_lat,lrc_lon ;
	FILE *filhisto ;


	// count hist is the monthly histogram for each alert

	int this_mday ;

	if (argc < 7) {
		cout << "Usage: modis_hist flist log_prefix modisflag ulc_lon ulc_lat lrc_lon lrc_lat " << endl ;
		exit(-1) ;
	}

	//command line arguments
	strcpy (flist, *++argv) ;
	strcpy (prefix, *++argv) ;
	modisflag = atoi (*++argv) ;
	ulc_lon = atof (*++argv) ;
	ulc_lat = atof (*++argv) ;
	lrc_lon = atof (*++argv) ;
	lrc_lat = atof (*++argv) ;


	modis_hdf *therm, *geo ;

	float *mnstdv ;


	// input file is an ascii file, each line has the radiance and the
	// geom file
	// input file contains the full path names of each file 

	// granule dimensions
	ns_modis = 1354 ;
	nl_modis = 2030 ;
	npix_modis = ns_modis * nl_modis ;

	//
	// hawaii
	//cent_lat = 19.4 ;
	//cent_lon = -155.25 ;
	// iceland
	//cent_lat = 64.64 ;
	//cent_lon = -17.528 ;


	FILE *ffile = fopen (flist, "r") ;


	sinu_1km *sinu  = new sinu_1km () ;
	alert *al = new alert() ;

	while (!feof (ffile)) {
	fscanf (ffile, "%s %s", infile, infile_geo) ; 
	tmpname = strstr (infile, "021KM")-3 ;

	therm = new modis_hdf (infile, modisflag) ;
	therm->get_date_period (infile, datearr) ;
	unxtime = converttime_unix (datearr[0], datearr[1], datearr[2], datearr[3], datearr[4], 0) ;
	thismonth = datearr[1] - 1 ;
	this_mday = mon_days[thismonth];
	geo = new modis_hdf (infile_geo,modisflag ) ;
	geo->get_date_period (infile_geo, datearr) ;
	geo->init_MOD03() ;

	if (geo->dayflag) {
		cout << "infile  : " << infile << endl << "DAYFILE" << endl ;
		delete therm ;
		delete geo ;
		continue ;
	}

	if (!geo->geom_status) {
		//cout << infile_geo << "UHOHOOH"<< endl ;
		cout << "problem with MOD03 FILE" << endl ;
	}

	filecount++ ;
	//if (filecount >10) break ;

	strcpy (fhisto, prefix) ;
	sprintf (fhisto, "%s_histo_%03d.txt", fhisto, this_mday) ;
	filhisto = fopen (fhisto, "a+") ;

	
	lat = geo->latarr ;
	lon = geo->lonarr ;
	sinu->set_badpix(geo->badpix) ;
	b21 = &therm->raddata_cal[0] ;
	b22 = &therm->raddata_cal[npix_modis] ;
	b32 = &therm->raddata_cal[3*npix_modis] ;

	//gridnum = sinu->get_gridnum(20., -157, &xval, yval) ;

	// get 4micron mean and stddev for each pixel in the modis image
	// if (this_mday !=1) break ;
	sinu->set_mday (this_mday) ;
	// zero out alerthist
	//
	for (i=0; i<npix_modis; i++) {
		latval = lat[i] ;
		lonval = lon[i] ;
		// check within the box
		if (latval > ulc_lat || latval < lrc_lat) continue ;
		if (lonval < ulc_lon || lonval > lrc_lon) continue ;
		b221val = b22[i] ;
		b32val = b32[i] ;
		if (b221val > 2.) b221val = b21[i] ;
		Tb221val = al->bb_radtotemp(4, b221val) ;
		Tb32val = al->bb_radtotemp(12, b32val) ;
		fprintf (filhisto,"%f\t%f\t%f\t%f\r\n", lonval, latval, Tb221val, Tb32val) ;
		fflush(filhisto) ;
	}

	



	delete therm ;
	delete geo ;
	} // read in the next pair of modis scenes

	delete al ;
	delete sinu ;
	fclose (ffile) ;
	fclose (filhisto) ;


}
