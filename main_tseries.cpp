#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "modis_hdf.h"


int main (int argc, char *argv[]) {

	char ifile[240], ifile_geo[240], ifil_list[240] ;
	char istring[420], *tmpname, daystr[40]  ;
	int i, mini, datearr[5], npix_modis ;
	modis_hdf *therm, *geom ;
	float lat0, lon0, latval, lonval, xdist, ydist, dist, dayfrac, day ; 
	float xdegm, ydegm, mindist ;
	float *b22, *b32, b22val, b32val ;
	ydegm = 110000. ;

	lon0 = -26.504 ;
	lat0 = 73.968 ;
	npix_modis = 1354 * 2030 ;
	xdegm = cos(lat0 * .017453) * ydegm ;

	strcpy (ifil_list, "filelist_iceland_terra.txt") ;
	FILE *flog = fopen ("tseries.txt", "w") ;
	FILE *fin = fopen (ifil_list,"r") ;
	while (!feof(fin)) {
		fscanf (fin, "%s %s", ifile, ifile_geo) ;
		tmpname = strstr (ifile,"021KM") - 3 ;
		geom = new modis_hdf (ifile_geo, 3) ;
		therm = new modis_hdf (ifile, 3) ;
		therm->get_date_period (ifile, datearr) ;
		geom->get_date_period (ifile_geo, datearr) ;
		geom->init_MOD03() ;

		mini = -1 ;
		mindist = 5000 ;

		for (i=0; i<1354L*2030; i++) {
			
			if (geom->badpix[i]==1) continue ;
			latval = geom->latarr[i] ;
			lonval = geom->lonarr[i] ;
			xdist = (lonval - lon0) * xdegm ;
			ydist = (latval - lat0) * ydegm ; 
			dist = sqrt (xdist * xdist + ydist * ydist) ;
			// check if within 5km
			if (dist < mindist) {
				mindist = dist;
				mini = i ;
			}
		}

		if (mindist < 5000) {
			b22 = &therm->raddata_cal[npix_modis] ;
			b32 = &therm->raddata_cal[2*npix_modis] ;
			b22val = b22[mini] ;
			b32val = b32[mini] ;
			strncpy (daystr, &tmpname[14], 3) ;
			dayfrac = datearr[3] / 24. + datearr[4]/1440. ;
			daystr[4]='\0' ;
			day = atoi(daystr) + dayfrac ;
			cout << "********************************"<< endl ;
			cout << ifile << " \t"<< day<< " \t" << b22val << " \t"<< b32val << endl ;
			fprintf (flog, "%s \t %f \t %f \t %f\r\n", ifile, day, b22val, b32val) ; 
			fflush(flog) ;
		}
		delete geom ;
		delete therm ;
	}
	fclose (flog) ;
	fclose (fin) ;

}


				
				
				



