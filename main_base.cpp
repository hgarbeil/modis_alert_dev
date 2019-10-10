#include "sinu_1km.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



//usage : makebase outfile month_num(0=Jan) 

int main (int argc, char *argv[]) {

	char outfile[420] ;
	float start_lat, start_lon, gspace, *lat, *lon ;
	float *mnstdv ;
	int   i, j, nx, ny, mnum ;
	int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;
	unsigned char *badpix ;

	strcpy (outfile, *++argv) ;
	mnum = atoi (*++argv) ;

	nx = 3800 ;
	ny = 3500 ;
	mnstdv = new float [nx * ny * 2] ;
	gspace = .008 ;
	start_lat = 64.64 + ny / 2 * gspace ;
	start_lon = -17.528 - nx / 2 * gspace ;

	lat = new float [nx * ny] ;
	lon = new float [nx * ny] ;
	badpix = new unsigned char [nx * ny] ;

	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) {
			lat[i*nx+j] = start_lat - i * gspace ;
			lon[i*nx+j] = start_lon + j * gspace ;
			badpix[i*nx+j] = 0 ;
		}
	}
		


	sinu_1km *sinu  = new sinu_1km() ;
	sinu->set_badpix(badpix) ;
	
	sinu->set_mday (mon_days[mnum]) ;
	sinu->get_mn_stdev_array (lat, lon, nx * ny, 2, mnstdv) ;

	FILE *fout = fopen (outfile, "w") ;
	fwrite (mnstdv, 4, 2 * nx * ny, fout) ;
	fclose (fout) ;

	delete [] lat ;
	delete [] lon ;
	delete [] mnstdv ;
}
