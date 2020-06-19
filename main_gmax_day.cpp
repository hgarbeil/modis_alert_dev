/***
* main_gma.cpp
* similar to the main_max program but instead of using the MOD021KM granule, it uses Eric's 4 banded radiance
* therefore we have a modis_hdf class 
* compares the incoming pixel 4 micron radiance to the global max value. This information currently resides 
* in the /local/worldbase/021km/full directory with files having a naming of rad221_highest_vvv_xx_yy_zz.bsq.
* vvv - is 1st day of the month (1, 32...
* xx - is modisflag
* yy - is row number
*zz - is column number
***/
#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "modis_hdf.h"
#include "modis_rad.h"
#include "modis_process.h"
#include "alert.h"
#include "testime.h"

#define DTOR 0.01745329252

extern int mon_days[] ;

int main (int argc, char *argv[]) {

	char webdir[420], webexc[240],  webexcc[240], flist[420], infile [420], *tmpname,  prefix[420], logfile[420] ;
	char infile_geo [420], daystr[4], outstr[240], tempstr[48]  ;
	char fwebfile[420] ;
	unsigned char *locdata, ucval, *almap ;
	
	float filecount =0, startlat, startlon, gspace ;
	float *lat, *lon, *b21, *b32, *b22, *b6, latval, lonval, xdist, ydist, dist ;
	int i, uval, thismonth, ns_modis, nl_modis, npix_modis, datearr[5], datearrfx[5], modisflag ;
	int nx, ny, npix_grid, gridnum, alerthist[7] ;
	bool firsttime, ntiflag, maxflag ;
	
	long unxtime ;

	// count hist is the monthly histogram for each alert
	int *monthhist = new int [12 * 4] ;
	for (i=0;i<48; i++) monthhist[i]  =0 ;

	int this_mday ;

	if (argc < 3) {
		cout << "Usage: modis_gmax flist  log_prefix modisflag" << endl ;
		exit(-1) ;
	}

    // get arguments from command line
	strcpy (flist, *++argv) ;
	strcpy (prefix, *++argv) ;
	modisflag = atoi (*++argv) ;
	strcpy (daystr, *++argv) ;
	strcpy (webdir,"/hotspot3/data2/maxalerts_day/2014/") ;
	strcpy (webexc,"/hotspot3/data2/alerts_exclusive_day/2014/") ;
	strcat (webdir, daystr) ;
	strcat (webexc, daystr) ;
	strcat (webdir, "/") ;
	strcat (webexc, "/") ;
	strcpy (webexcc, webexc) ;

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
	// alert files for alert website viewing
	FILE *fweb ;
	FILE *fwebex_nti, *fwebex_max ;

                    // write headings for ASCII column files
	fprintf (flog, "FName\tJDay\tUnixtime\t#NTI\tDEL35\tDEL40\tDEL45\r\n") ;
	fprintf (fstd_log, "FName\tUnixtime\tLon\tLat\tSamp\tLine\tNTI\tTDiff\tb221\tb22MAX\tFRAC\tb32max\tSolZenith\r\n") ;
	fprintf (fnti_log, "FName\tUnixtime\tLon\tLat\tSamp\tLine\tNTI\tDELTb\tB22divB22MAX\r\n") ;
	

	modis_hdf *geo ;
	modis_rad *therm ;
	

	float *max2232 ;
	vector <int> al_nti_inds ;
	vector <float> al_nti_vals ;
	vector <int> al_max_inds ;
	vector <float> al_max_vals ;
	alert *al = new alert() ;


	// input file is an ascii file, each line has the radiance and the
	// geom file
	// input file contains the full path names of each file 
                    // size of output tile (800m pixels)
	//nx = 3800 ;
	//ny = 3500 ;
	ns_modis = 1354 ;
	nl_modis = 2030 ;

	npix_modis = ns_modis * nl_modis ;

    // alloc arrays for input global mn, stdv and output arrays
	max2232 = new float [npix_modis *2] ;
	almap = new unsigned char [npix_modis] ;

	FILE *ffile = fopen (flist, "r") ;

	
	float *al_nti = new float [1354L * 2030] ;
	float *al_max = new float [1354L * 2030] ;
	float tdiff ;

	sinu_1km *sinu  = new sinu_1km () ;
	// DEBUG - write out basegrid
	// FILE *fout0 = fopen ("basegrid", "w") ;



	while (!feof (ffile)) {
	ntiflag = false ;
	maxflag = false ;
	fscanf (ffile, "%s %s", infile, infile_geo) ; 
	tmpname = strstr (infile, "021KM.A201")-3 ;
	fweb =NULL ;
	for (i=0; i<npix_modis; i++) almap[i]=0;

	geo = new modis_hdf (infile_geo,modisflag ) ;
	geo->get_date_period (infile_geo, datearr) ;
	geo->get_date_period_fix (infile_geo, datearrfx) ;
	geo->init_MOD03() ;
	// for daytime - check for glint angle and mark as bad pixels

	//exit(1) ;

	//unxtime = converttime_unix (datearr[0], datearrfx[1], datearrfx[2], datearrfx[3], datearrfx[4], 0) ;
	unxtime = converttime_unix (datearr[0], datearrfx[1], datearrfx[2], datearrfx[3], datearrfx[4], 0) ;
	thismonth = datearr[1] - 1 ;
	this_mday = mon_days[thismonth];
    // use the section below when the global dataset is not complete
	/*
	if (this_mday != 1 && this_mday <62 )  {
		delete therm ;
		continue ;
	}
	*/
	therm = new modis_rad (infile) ;
	therm->readfile() ;
	al->set_max (therm->b21max, therm->b22max) ;

	if (!geo->dayflag) {
		cout << "infile  : " << infile << endl << "NIGHTFILE" << endl ;
		delete therm ;
		delete geo ;
		continue ;
	}

	if (!geo->geom_status) {
		//cout << infile_geo << "UHOHOOH"<< endl ;
		cout << "problem with MOD03 FILE" << endl ;
	}
	// calc glint, mark glint < 12 as badpix
	geo->calc_glint() ;
	filecount++ ;
	//if (filecount >10) break ;

	//FILE *fglin = fopen ("/home/harold/workdir/modis_alert_dev/badpix", "w") ;
	//fwrite (geo->badpix, 1, npix_modis,fglin) ;
	//fclose (fglin) ;
	//exit(0) ;
	
	lat = geo->latarr ;
	lon = geo->lonarr ;
	al->set_badpix (geo->badpix) ;
	sinu->set_badpix(geo->badpix) ;
	b21 = therm->b21 ;
	b22 = therm->b22 ;
	b32 = therm->b32 ;
	b6 = therm->b6 ;


	// get 4micron mean and stddev for each pixel in the modis image
	sinu->set_mday (this_mday) ;
	//sinu->get_mn_stdev_array (lat, lon, 1354L*2030, modisflag, mnstdv) ;
	
    // load up the max2232 arrays , bsq with 22 in band 1 and 32 in band2
    sinu->get_max_array (lat, lon, 1354L*2030, modisflag, max2232) ;
    // DEBUG - output b221 mean and stdv
	
                   
	al_nti_inds.clear() ;
	al_nti_vals.clear() ;
	al_max_inds.clear() ;
	al_max_vals.clear() ;
	// zero out alerthist
                    // zero out the alert histogram
	memset (alerthist, 0, sizeof (alerthist)) ;

        
        
	// 0 is nti hits
	//int num_std_hits = al->calc_nstdv (b21, b22, mnstdv, al_max, al_max_inds, al_max_vals) ;
    int num_max_hits = al->calc_max_day_dev (b6, b21, b22, b32, max2232, al_max, al_max_inds, al_max_vals) ;
    //int num_max_hits = al->calc_max_dev (b21, b22, b32, max2232, al_max, al_max_inds, al_max_vals) ;
	int num_hits = al->calc_nti_day (b6, b21, b22, b32, al_nti, al_nti_inds, al_nti_vals) ; 
    //int num_max_hits = al->calc_max (b21, b22, b32, max2232, al_max, al_max_inds, al_max_vals) ;
	if (num_hits==0 && num_max_hits==0) {
		sprintf (outstr, "%s %d\t %ld\t %d\t %d\t %d\t %d \r\n", 
			tmpname, datearr[2], unxtime, alerthist[0], alerthist[1], alerthist[2], 
			alerthist[3]) ;
		fputs (outstr, flog) ;
		fflush (flog) ;
		delete therm ;
		delete geo ;
		continue ;
	}


	// now fill the layers grids 
	int index,  iline, isamp ;
	//cout << "ALERTS : " << al_nti_inds.size()  <<  endl ;
	//cout << "3std ALERTS : " << al_max_inds.size()  <<  endl ;
	//alerthist[0] = al_nti_inds.size() ;
	//alerthist[1] = al_max_inds.size() ;
	float b221, b32val ;
	int count = 0 ;


	// exclusive files
			strcpy (fwebfile, webdir) ; 
			strncpy (tempstr, tmpname, 22) ;
			tempstr[22]='\0';
			strcat (tempstr, "_web.txt") ;
			strcat (fwebfile, tempstr) ;
			fweb = fopen (fwebfile, "w") ;

			strcpy (fwebfile, webexc) ;
			strncpy (tempstr, tmpname, 22) ;
			tempstr[22]='\0';
			strcat (tempstr, "_ntix.txt") ;
			strcat (fwebfile, tempstr) ;
			fwebex_nti = fopen (fwebfile, "w") ;

			strcpy (fwebfile, webexcc) ;
			strncpy (tempstr, tmpname, 22) ;
			tempstr[22]='\0';
			strcat (tempstr, "_maxx.txt") ;
			strcat (fwebfile, tempstr) ;
			fwebex_max = fopen (fwebfile, "w") ;

			fprintf (fweb, "%04d-%02d-%02dT%02d:%02d:00\r\n",datearrfx[0],datearrfx[1],datearrfx[2],datearrfx[3],datearrfx[4]); 
			fprintf (fwebex_nti, "%04d-%02d-%02dT%02d:%02d:00\r\n",datearrfx[0],datearrfx[1],datearrfx[2],datearrfx[3],datearrfx[4]); 
			fprintf (fwebex_max, "%04d-%02d-%02dT%02d:%02d:00\r\n",datearrfx[0],datearrfx[1],datearrfx[2],datearrfx[3],datearrfx[4]); 


	for (i=0; i<al_nti_inds.size() ; i++) {
		index = al_nti_inds.at(i) ;
		iline = index / ns_modis ;
		isamp = index - (iline * ns_modis) ;
		latval = lat[index] ;
		lonval = lon[index] ;
		// if an nti hit mark with a 1
		almap[index] = 1 ;
		float fracval = b22[index] / max2232[index] ;
		fprintf (fnti_log, "%s %ld %f %f %d %d %f %f %f\r\n", tmpname, unxtime,
			lonval, latval, isamp, iline, al_nti[index], al_max[index], fracval) ;
		fflush (fnti_log);
		alerthist[0]++ ;
		ntiflag = true ;
		
	}
	for (i=0; i<al_max_inds.size() ; i++) {
		index = al_max_inds.at(i) ;
		tdiff = al_max[index] ;

		// 35 deg TDiff
		if (al_max_vals.at(i) < .999 || tdiff < 15) continue ;
		// if an nti hit, both alerts are marked with a 3
		// if only a max hit, then its marked with a 2
		if (almap[index]==1) 
			almap[index]= 3 ;
		else 
			almap[index] = 2 ;


		if (max2232[index] <= .0001) {
			int hg =1 ;
			hg *=1 ;
		}
		iline = index / ns_modis ;
		isamp = index - (iline * ns_modis) ;
		latval = lat[index] ;
		lonval = lon[index] ;
		if (al_max_vals.at(i) < 1.05) continue ; 
		
		if (tdiff > 15.) {
			alerthist[1]++ ;
		}
		if (tdiff > 20.) {
			alerthist[2]++ ;
		}
		if (tdiff>25.) {
			alerthist[3]++ ;
		}

		b221 = b22[index] ;
		b32val = b32[index] ;
		if (b221 > 2.0) b221 = b21[index] ;
		if (b221 < -9) b221 = b21[index] ;
		//if (b221 < 0.0001) b221 = 0. ;
		fprintf (fstd_log, "%s %ld %f %f %d %d %f %f %f %f %f %f %f\r\n", tmpname, unxtime,
			lonval, latval, isamp, iline, al_nti[index], al_max[index],  b221, max2232[index], b221/max2232[index], max2232[index+npix_modis], geo->solzen[index]/100.) ; 
		fflush (fstd_log);

		fprintf (fweb, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %9.4f %7.4f %7.4f %7.4f %4d %4d %8.4f %8.4f\r\n",
			b21[index], b22[index], b6[index], b32[index], b32[index], latval, lonval,  
			geo->solzen[index]/100.*DTOR, geo->senszen[index]/100.*DTOR, geo->sensaz[index]/100.*DTOR, iline, isamp, al_max[index], b221/max2232[index]) ;
		fflush (fweb) ;
		maxflag = true ;
	}
			
	sprintf (outstr, "%s %d\t %ld\t %d\t %d\t %d\t %d\r\n", tmpname, datearr[2], unxtime, alerthist[0], alerthist[1], alerthist[2], 
		alerthist[3]) ;
	fputs (outstr, flog) ;
	fflush (flog) ;

	// now write out exclusive alerts 
	for (i=0; i<npix_modis; i++) {
		b221 = b22[i] ;
		b32val = b32[i] ;
		iline = i/ ns_modis ;
		isamp = i- (iline * ns_modis) ;
		if (b221 > 2.0) b221 = b21[index] ;
		//if (b221 < 0.) b221 = b21[index] ;
		if (b221 < 0.001 && b221 > -1) continue ;
		if (b221 < -9.9) b221 = b21[index] ;
	
		if (almap[i]==1) { 
			fprintf (fwebex_nti, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %9.4f %7.4f %7.4f %7.4f %4d %4d %7.4f %8.4f %8.4f\r\n",
			b21[i], b22[i], b6[i], b32[i], b32[i], lat[i], lon[i],  
			geo->solzen[i]/100.*DTOR, geo->senszen[i]/100.*DTOR, geo->sensaz[i]/100.*DTOR, iline, isamp, al_nti[i], al_max[i], b221/max2232[i]) ;

			fflush (fwebex_nti) ;
		}
		if (almap[i]==2) { 
			fprintf (fwebex_max, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %9.4f %7.4f %7.4f %7.4f %4d %4d %7.4f %8.4f %8.4f\r\n",
			b21[i], b22[i], b6[i], b32[i], b32[i], lat[i], lon[i],  
			geo->solzen[index]/100.*DTOR, geo->senszen[i]/100.*DTOR, geo->sensaz[i]/100.*DTOR, iline, isamp, al_nti[i], al_max[i], b221/max2232[i]) ;

			fflush (fwebex_max) ;
		}
	}


	if (fweb)
	{
		fclose (fweb) ;
		fclose (fwebex_max) ;
		fclose (fwebex_nti) ;
	}



	cout << "Number of old hot spots is " << alerthist[0] << endl ;
	cout << "Number of new hot spots is " << alerthist[1] << endl ;
	delete therm ;
	delete geo ;

	} // read in the next pair of modis scenes

	fclose (fstd_log) ;
	fclose (fnti_log) ;
	fclose (flog) ;
	//fclose (fweb) ;
        
        
	delete [] al_nti ;
	delete [] al_max ;
	delete [] max2232 ;

}
