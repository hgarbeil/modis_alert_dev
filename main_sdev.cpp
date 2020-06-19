#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "modis_hdf.h"
#include "modis_process.h"
#include "alert.h"
#include "testime.h"

extern int mon_days[] ;

int main (int argc, char *argv[]) {

	char flist[420], infile [420], *tmpname,  prefix[420], logfile[420] ;
	char flocfile [420] ;
	char monlog[420], flayers[420] ;
	char infile_geo [420], outstr[240]  ;
	unsigned char *locdata, ucval ;
	
	float filecount =0, startlat, startlon, gspace ;
	float *lat, *lon, *b21, *b32, *b22, latval, lonval, xdist, ydist, dist ;
	int i, uval, thismonth, ns_modis, nl_modis, npix_modis, datearr[5], datearrfx[5], modisflag ;
	int nx, ny, npix_grid, gridnum, alerthist[7] ;
	long unxtime ;
	float cent_lat, cent_lon ;

	// count hist is the monthly histogram for each alert
	int *monthhist = new int [12 * 4] ;
	for (i=0;i<48; i++) monthhist[i]  =0 ;

	int this_mday ;

	if (argc < 3) {
		cout << "Usage: modis_sdev flist  log_prefix modisflag c_lon c_lat" << endl ;
		exit(-1) ;
	}

                    // get arguments from command line
	strcpy (flist, *++argv) ;
	strcpy (prefix, *++argv) ;
	modisflag = atoi (*++argv) ;
	cent_lon = atof (*++argv) ;
	cent_lat = atof (*++argv) ;

                    // input files 
	strcpy (logfile, prefix) ;
	strcat (logfile, "_log.txt") ;
	FILE *flog = fopen (logfile, "w") ;
	strcpy (logfile, prefix) ;
	strcat (logfile, "_nti.txt") ;
	FILE *fnti_log = fopen (logfile, "w") ;
	strcpy (logfile, prefix) ;
	strcat (logfile, "_std.txt") ;
	FILE *fstd_log = fopen (logfile, "w") ;
	strcpy (monlog, prefix) ;
	strcat (monlog, "_month.txt") ;
	FILE *fmon_log = fopen (monlog, "w") ;
	fprintf (flog, "FName\tJDay\tUnixtime\t#NTI\t#3_SDEV\t#4_SDEV\t#5_SDEV\r\n") ;
	fprintf (fstd_log, "FName\tUnixtime\tLon\tLat\tSamp\tLine\tNTI\t#Stdvs\tb221\tb22mean\tb22stdv\tdist\tSolZenith\r\n") ;
	fprintf (fnti_log, "FName\tUnixtime\tLon\tLat\tSamp\tLine\tNTI\t#Stdvs\r\n") ;
	fprintf (fmon_log, "Month\t#NTI\t#3_SDEV\t#4_SDEV\t#5_SDEV\r\n") ;
	strcpy (flayers, prefix) ;
	strcat (flayers, "_layers") ;
	strcpy (flocfile, prefix) ;
	strcat (flocfile, "_location") ;

	modis_hdf *therm, *geo ;

	float *mnstdv ;
	vector <int> al_nti_inds ;
	vector <float> al_nti_vals ;
	vector <int> al_std_inds ;
	vector <float> al_std_vals ;
	alert *al = new alert() ;


	// input file is an ascii file, each line has the radiance and the
	// geom file
	// input file contains the full path names of each file 
                    // size of output tile (800m pixels)
	//nx = 3800 ;
	//ny = 3500 ;
	nx = 2048 ;
	ny = 2048 ;
	ns_modis = 1354 ;
	nl_modis = 2030 ;

	// hawaii
	//cent_lat = 19.4 ;
	//cent_lon = -155.25 ;
	// iceland
	//cent_lat = 64.64 ;
	//cent_lon = -17.528 ;
                  // hawaii
	//cent_lat = 19.4 ;
	//cent_lon = -155.25 ;
	// iceland
	//cent_lat = 64.64 ;
	//cent_lon = -17.528 ;
	// newguinea
	//cent_lat = 0. ;
	//cent_lon = 138. ;21
	npix_modis = ns_modis * nl_modis ;

                    // alloc arrays for input global mn, stdv and output arrays
	mnstdv = new float [npix_modis *2] ;
	locdata = new unsigned char [nx * ny * 3] ;

	

	FILE *ffile = fopen (flist, "r") ;

	//strcpy (infile, "/hotspot3/data2/modis/hawaii/terra/MOD021KM.A2014003.0845.061.2017307163855.hdf") ;
	//strcpy (infile_geo, "/hotspot3/data2/modis/hawaii/terra/MOD03.A2014003.0845.061.2017307163343.hdf") ;


	
	npix_grid = nx * ny ;
	gspace = .008 ;
	startlat =  cent_lat + ny/2. * gspace ;
	startlon = cent_lon - nx/2 * gspace ;
	float *basegrid = new float [2*nx * ny] ;
	unsigned short *layers = new unsigned short  [6 * nx * ny] ;
	float *al_nti = new float [1354L * 2030] ;
	float *al_std = new float [1354L * 2030] ;

	sinu_1km *sinu  = new sinu_1km () ;
	sinu->makegrid (basegrid, startlat, startlon, gspace, nx, ny) ;
	// DEBUG - write out basegrid
	// FILE *fout0 = fopen ("basegrid", "w") ;

	for (i=0; i<npix_grid; i++) {
		layers[i] = basegrid[i*2]*100 ;
		uval = 1 * layers[i] ;
		if (uval < 0) uval = 0 ;
		if (uval > 100) uval = 100 ;
		ucval = (unsigned char) uval ;
		layers[npix_grid+i] = basegrid[i*2+1]*1000 ;
		locdata[i] = ucval ;
		locdata[npix_grid+i] = ucval ;
		locdata[2*npix_grid+i] = ucval ;
		layers[2*npix_grid+i] = 0 ; 
		layers[3*npix_grid+i] = 0; 
		layers[4*npix_grid+i] = 0; 
		layers[5*npix_grid+i] = 0; 
	}	


	while (!feof (ffile)) {
	fscanf (ffile, "%s %s", infile, infile_geo) ; 
	tmpname = strstr (infile, "021KM")-3 ;

	therm = new modis_hdf (infile, modisflag) ;
	therm->get_date_period (infile, datearr) ;
	therm->get_date_period_fix (infile, datearrfx) ;
	unxtime = converttime_unix (datearrfx[0], datearrfx[1], datearrfx[2], datearrfx[3], datearrfx[4], 0) ;
	al->set_max (therm->b21max, therm->b22max) ;
	thismonth = datearr[1] - 1 ;
	this_mday = mon_days[thismonth];
	/*
	if (this_mday != 1 && this_mday <62 )  {
		delete therm ;
		continue ;
	}
	*/
	if (this_mday < 184) {
		delete therm ;
		continue ;
	}
	if (this_mday > 185) {
		delete therm;
		break ;
	}
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
	al->set_badpix (geo->badpix) ;
	sinu->set_badpix(geo->badpix) ;
	b21 = &therm->raddata_cal[0] ;
	b22 = &therm->raddata_cal[npix_modis] ;
	b32 = &therm->raddata_cal[3*npix_modis] ;

	//gridnum = sinu->get_gridnum(20., -157, &xval, yval) ;

	// get 4micron mean and stddev for each pixel in the modis image
	sinu->set_mday (this_mday) ;
	sinu->get_mn_stdev_array (lat, lon, 1354L*2030, modisflag, mnstdv) ;
	// DEBUG - output b221 mean and stdv
	/*
	FILE *fout2 = fopen ("mnstdv", "w") ;
	fwrite (mnstdv, 8, 1354L * 2030, fout2) ;
	fclose (fout2) ;
	*/
	// DEBUG - output b22 
	/*
	FILE *fout3 = fopen ("b221", "w") ;
	fwrite (b21, 4, 1354L * 2030, fout3) ;
	fwrite (b22, 4, 1354L * 2030, fout3) ;
	fclose (fout3) ;
	*/
                    //void getNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *outarr, bool openFlag) ;
                    //sinu->getNTIValues (thismonth, 0, lat, lon, 1354*2030,b22_5, nti_5, 0) ;
               
	al_nti_inds.clear() ;
	al_nti_vals.clear() ;
	al_std_inds.clear() ;
	al_std_vals.clear() ;
	// zero out alerthist
	memset (alerthist, 0, sizeof (alerthist)) ;

        
        
	int num_hits = al->calc_nti (b21, b22, b32, al_nti, al_nti_inds, al_nti_vals) ; 
	// 0 is nti hits
	int num_std_hits = al->calc_nstdv (b21, b22, mnstdv, al_std, al_std_inds, al_std_vals) ;
	if (num_hits==0 && num_std_hits==0) {
		sprintf (outstr, "%s %d\t %ld\t %d\t %d\t %d\t %d \r\n", tmpname, datearr[2], unxtime, alerthist[0], alerthist[1], alerthist[2], 
			alerthist[3]) ;
		fputs (outstr, flog) ;
		fflush (flog) ;
		delete therm ;
		delete geo ;
		continue ;
	}


	// now fill the layers grids 
	int index, grid_x, grid_y, iline, isamp ;
	//cout << "ALERTS : " << al_nti_inds.size()  <<  endl ;
	//cout << "3std ALERTS : " << al_std_inds.size()  <<  endl ;
	//alerthist[0] = al_nti_inds.size() ;
	//alerthist[1] = al_std_inds.size() ;
	float b221 ;
	int count = 0 ;
	for (i=0; i<al_nti_inds.size() ; i++) {
		index = al_nti_inds.at(i) ;
		iline = index / ns_modis ;
		isamp = index - (iline * ns_modis) ;
		latval = lat[index] ;
		lonval = lon[index] ;
		grid_x = int((lonval - startlon) /gspace +0.5) ; 
		grid_y = int((startlat - latval) /gspace +0.5) ; 
		layers[2*npix_grid + grid_y*nx+grid_x]++ ;
		
		fprintf (fnti_log, "%s %ld %f %f %d %d %f %f\r\n", tmpname, unxtime,
			lonval, latval, isamp, iline, al_nti[index], al_std[index]) ;
		fflush (fnti_log);
		alerthist[0]++ ;
		
	}
	for (i=0; i<al_std_inds.size() ; i++) {
		index = al_std_inds.at(i) ;
		iline = index / ns_modis ;
		isamp = index - (iline * ns_modis) ;
		latval = lat[index] ;
		lonval = lon[index] ;
		grid_x = int((lonval - startlon) /gspace +0.5) ; 
		grid_y = int((startlat - latval) /gspace +0.5) ; 
		if (grid_x < 0 || grid_y < 0) continue ;
		if (grid_x >= nx || grid_y >= ny) continue ;
		if (al_std_vals.at(i) > 3.0) {
			layers[3 * npix_grid + grid_y*nx+grid_x]++ ;
			// 1 is 2 sdev
			monthhist[thismonth*4+1]++ ;
			locdata[grid_y*nx+grid_x] = 255 ;
			locdata[npix_grid + grid_y*nx+grid_x] = 0 ;
			locdata[2*npix_grid + grid_y*nx+grid_x] = 0 ;
			alerthist[1]++ ;
		}
		if (al_std_vals.at(i) > 4) {
			layers[4 * npix_grid + grid_y*nx+grid_x]++ ;
			alerthist[2]++ ;
			// 2 is 2.6 sdev
			monthhist[thismonth*4+2]++ ;
			locdata[grid_y*nx+grid_x] = 0 ;
			locdata[npix_grid + grid_y*nx+grid_x] = 255 ;
			locdata[2*npix_grid + grid_y*nx+grid_x] = 0 ;
		}
		if (al_std_vals.at(i) > 5.) {
			layers[5 * npix_grid + grid_y*nx+grid_x]++ ;
			alerthist[3]++ ;
			// 3 is 3 sdev
			monthhist[thismonth*4+3]++ ;
			locdata[grid_y*nx+grid_x] = 0 ;
			locdata[npix_grid + grid_y*nx+grid_x] = 0 ;
			locdata[2 * npix_grid + grid_y*nx+grid_x] = 255 ;
		}

		cout << i << " : " << index << "   " << grid_x<< " , " << grid_y << endl ;
		b221 = b22[index] ;
		ydist = latval - cent_lat ;
		xdist = lonval - cent_lon ;
		dist = 100. * (sqrt (xdist * xdist + ydist * ydist)) ;
		if (b221 > 1.9) b221 = b21[index] ;
		fprintf (fstd_log, "%s %ld %f %f %d %d %f %f %f %f %f %f %f\r\n", tmpname, unxtime,
			lonval, latval, isamp, iline, al_nti[index], al_std[index],  b221, mnstdv[index*2], mnstdv[index*2+1], dist,
				geo->solzen[index]/100.) ; 
		fflush (fstd_log);
	}
			
	sprintf (outstr, "%s %d\t %ld\t %d\t %d\t %d\t %d\r\n", tmpname, datearr[2], unxtime, alerthist[0], alerthist[1], alerthist[2], 
		alerthist[3]) ;
	fputs (outstr, flog) ;
	fflush (flog) ;


	cout << "Number of old hot spots is " << alerthist[0] << endl ;
	cout << "Number of new hot spots is " << alerthist[1] << endl ;
	delete therm ;
	delete geo ;

	} // read in the next pair of modis scenes

	for (i=0; i<12; i++) {
		sprintf (outstr, "%d\t%d\t%d\t%d\t%d\r\n",  
			i+1, monthhist[i*4], monthhist[i*4+1], monthhist[i*4+2],
			monthhist[i*4+3]) ;
		fputs (outstr, fmon_log) ;
	}
	fclose (fmon_log) ;
	fclose (fstd_log) ;
	fclose (fnti_log) ;
	fclose (flog) ;
	FILE *fout = fopen (flayers, "w") ;
	char ohdr[240] ;
	strcpy (ohdr, flayers) ;
	strcat (ohdr, ".hdr") ;
	al->write_envi_header (ohdr, nx, ny, startlat, startlon, gspace) ;
	fwrite (layers, 2, 6 * npix_grid, fout) ;
	fclose (fout) ;
	strcpy (ohdr, flocfile) ;
	strcat (ohdr, ".hdr") ;
	al->write_envi_header (ohdr, nx, ny, startlat, startlon, gspace, 1) ;
	fout = fopen (flocfile, "w") ;
	fwrite (locdata, 1, 3 * npix_grid, fout) ;
	fclose (fout) ;

        
        
        
	delete [] al_nti ;
	delete [] al_std ;
	delete [] layers ;
	delete [] basegrid ;
	delete [] mnstdv ;
	delete [] locdata ;
	delete [] monthhist ;

}
