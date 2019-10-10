#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "modis_hdf.h"
#include "modis_process.h"
#include "alert.h"
#include "testime.h"

int main (int argc, char *argv[]) {

	char flist[420], infile [420], *tmpname,  prefix[420] ;
	char monlog[420], flayers[420], fcount[420],fbase[420] ;
	char infile_geo [420], outstr[240]  ;
	unsigned char *countdata; 
	float *histodata  ; 
	
	float filecount =0, startlat, startlon, gspace ;
	float *lat, *lon, *b21, *b31, *b22, latval, lonval, xdist, ydist, dist ;
	int i, countnum, thismonth, ns_modis, nl_modis, npix_modis, datearr[5], modisflag ;
	int lineloc, samploc, newloc ;
	int nx, ny, npix_grid, gridnum ;
	long unxtime ;
	float cent_lat, cent_lon ;

	int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;
	// count hist is the monthly histogram for each alert

	int this_mday ;

	if (argc < 3) {
		cout << "Usage: modis_sdev flist  log_prefix modisflag c_lon c_lat" << endl ;
		exit(-1) ;
	}


	strcpy (flist, *++argv) ;
	strcpy (prefix, *++argv) ;
	modisflag = atoi (*++argv) ;
	cent_lon = atof (*++argv) ;
	cent_lat = atof (*++argv) ;

	strcpy (flayers, prefix) ;
	strcpy (fcount, prefix) ;
	strcat (flayers, "_histo") ;
	strcat (fcount, "_count") ;

	modis_hdf *therm, *geo ;

	float *mnstdv ;


	// input file is an ascii file, each line has the radiance and the
	// geom file
	// input file contains the full path names of each file 

	// nx, ny will be size of histo file and basemap
	nx = 2400 ;
	ny = 2400 ;
	ns_modis = 1354 ;
	nl_modis = 2030 ;
	// hawaii
	//cent_lat = 19.4 ;
	//cent_lon = -155.25 ;
	// iceland
	//cent_lat = 64.64 ;
	//cent_lon = -17.528 ;
	npix_modis = ns_modis * nl_modis ;

	histodata = new float  [nx * ny * 90] ;
	countdata = new unsigned char [nx * ny] ;



	

	FILE *ffile = fopen (flist, "r") ;

	//strcpy (infile, "/hotspot3/data2/modis/hawaii/terra/MOD021KM.A2014003.0845.061.2017307163855.hdf") ;
	//strcpy (infile_geo, "/hotspot3/data2/modis/hawaii/terra/MOD03.A2014003.0845.061.2017307163343.hdf") ;


	
	npix_grid = nx * ny ;
	gspace = .008 ;
	startlat =  cent_lat + ny/2. * gspace ;
	startlon = cent_lon - nx/2 * gspace ;
	float *basegrid = new float [4*nx * ny] ;

	sinu_1km *sinu  = new sinu_1km () ;
	sinu->makegrid (basegrid, startlat, startlon, gspace, nx, ny) ;
	sinu->makegrid_32 (&basegrid[2*npix_grid], startlat, startlon, gspace, nx, ny) ;

	// initialize arrays to zero
	for (i=0; i<npix_grid*90L; i++)histodata[i] = 0 ; 
	for (i=0; i<npix_grid; i++) {
		countdata[i] = 0 ;
	}	


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

	
	lat = geo->latarr ;
	lon = geo->lonarr ;
	sinu->set_badpix(geo->badpix) ;
	b21 = &therm->raddata_cal[0] ;
	b22 = &therm->raddata_cal[npix_modis] ;
	b31 = &therm->raddata_cal[2*npix_modis] ;

	//gridnum = sinu->get_gridnum(20., -157, &xval, yval) ;

	// get 4micron mean and stddev for each pixel in the modis image
	if (this_mday !=1) break ;
	sinu->set_mday (this_mday) ;
	// zero out alerthist
	//
	for (i=0; i<npix_modis; i++) {
		lineloc = int((startlat - lat[i] )/gspace + 0.5) ;
		samploc = int((lon[i] - startlon) / gspace + 0.5) ;
		if (samploc < 0 || samploc > nx) continue ;
		if (lineloc < 0 || lineloc > ny) continue ;
		newloc = lineloc * nx + samploc ;
		countnum = countdata[newloc] ;
		histodata[countnum * npix_grid + newloc] = b22[i] ;
		countdata[newloc]++ ;
	}

	



	delete therm ;
	delete geo ;
	} // read in the next pair of modis scenes

	FILE *fout = fopen (flayers, "w") ;
	char ohdr[240] ;
	strcpy (ohdr, flayers) ;
	strcat (ohdr, ".hdr") ;
	//al->write_envi_header (ohdr, nx, ny, startlat, startlon, gspace) ;
	fwrite (histodata, 4, 90 * npix_grid, fout) ;
	fclose (fout) ;

	strcpy (ohdr, fcount) ;
	strcat (ohdr, ".hdr") ;
	//al->write_envi_header (ohdr, nx, ny, startlat, startlon, gspace, 1) ;
	fout = fopen (fcount, "w") ;
	fwrite (countdata, 1, 1 * npix_grid, fout) ;
	fclose (fout) ;


	delete [] basegrid ;
	delete [] countdata ;
	delete [] histodata ;

}
