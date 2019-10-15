 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sinu_1km.cpp
 * Purpose : Various functions for the sinu_1km class. This will relate the latitude and longitude 
 * based to an xy cartesian system based upon the sinusoidal projection. 
 * Author: harold garbeil
 * 
 * Created on August 3, 2018, 9:29 AM
 */

#include "sinu_1km.h"
#include "surftemp.h"

sinu_1km::sinu_1km() {
    Radius = 6371007.181 ;
    cent_merid = 0. ;
    mDegToRad = 0.01745329 ;
	meters_per_degree = 2. * M_PI * Radius / 360. ;
    wrld_starty =  meters_per_degree * 90. ; // +90 * metersperdeg
    wrld_startx = -2. * wrld_starty ; //-180 * metersperdeg
    gspace = 10. * meters_per_degree / 1200. ;
	gspace10 = 10. * meters_per_degree ;
	mday = 1 ;
	badpix = 0L ;
}

void sinu_1km::set_badpix (unsigned char *b) {
	badpix =b ;
}

void sinu_1km::set_mday (int md) {
	mday = md ;
}

// grid number divides latlon grid into 10 degree blocks
int sinu_1km::get_gridnum (float lat, float lon, int *xloc, int *yloc) {
	int iy, ix, igrid ;
	float x=0. , y=0.  ;
	float startx, starty, xval, yval ;
	latlon_to_xy (lat, lon, &x, &y)  ;

	
	xval = (x - wrld_startx) / gspace10 ;
	yval = (wrld_starty - y) / gspace10 ;
	ix = int (xval) ;
	iy = int (yval) ;
	startx = wrld_startx + gspace10 * ix ;
	starty = wrld_starty - gspace10 * iy ;
	*xloc = int((x - startx) / gspace) ;
	*yloc = int((starty - y) / gspace) ;


	if (ix < 0) ix = 0 ;
	if (iy < 0) iy = 0 ;
	if (iy > 17) iy = 17 ;
	if (ix > 35) ix = 35 ;
	igrid = iy * 36 + ix ;
	if (igrid > 500){
		cout<<"igrid is " << igrid << endl ;
		cout << "Lat : "<< lat<< endl ;
		cout << "Lon : "<< lon<< endl ;
		return (-1) ;
	}
		
	return  (iy * 36 + ix) ;

}

/* based on lat lon start and gspace, and nx, ny, fill the grid with
 * the mean and standard deviation band 21/22 value
 */
void sinu_1km::makegrid_32 (float *basegrid, float startlat, float startlon, float gspace, int nx, int ny) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [nx * ny * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<nx * ny * 3; i++) indii[i]=-1 ;

	for (i=0; i<ny; i++) {
		latval = startlat - i * gspace ;
		for (j=0; j<nx ; j++) {
			lonval = startlon + j * gspace ;
			gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

			if (firstgrid) {
				listgrids[ngrids++]=gridnum ;
				firstgrid=0 ;
				indii [i * nx * 3 + j * 3] = 0 ;
			}
			else {
				foundgrid=0 ;
				// not a new grid
				for (igrid=0; igrid<ngrids; igrid++) {
					if (gridnum == listgrids[igrid]) {
						indii [i * nx * 3 + j * 3] = igrid ;
						foundgrid=1 ;
						break ;
					}
				}
				// if a new grid - unique number
				if (foundgrid==0) {
					indii [i * nx * 3 + j * 3] = ngrids ;
					listgrids[ngrids++] = gridnum ;
				}
			}
					
			//indii[i*nx*3+j*3]=gridnum ;
			indii[i*nx*3+j*3+1]=xloc ;
			indii[i*nx*3+j*3+2]=yloc ;
		}
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L * ngrids] ;
	float *indata_stdv = new float [1200L * 1200L * ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		// note that these are bsq files, mean, then stdev
		sprintf (infile, "/local/worldbase/021km/orig_rad32_low_%03d_03_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_03_%02d_%02d.bsq\0",
			mday, irow, icol) ;
		fin = fopen (infile, "r") ;
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdv[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < nx * ny ; isamp++) {
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		basegrid[isamp*2] = indata [inloc] ;
		basegrid[isamp*2+1] = indata_stdv [inloc] ;
		//mnstdev[isamp*2+1] = indata_stdv [inloc] ;
	}
/*
	FILE *fout =fopen ("testfile", "w") ;
	fwrite (basegrid, 4,  2*nx*ny, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdv ;
	delete [] indii ;



}
	
/* given a lat and a lon array, this function returns the 4 micron (221) mean and stand deviation 
/* based on lat lon start and gspace, and nx, ny, fill the grid with
 * the mean and standard deviation band 21/22 value
 */
void sinu_1km::makegrid (float *basegrid, float startlat, float startlon, float gspace, int nx, int ny) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [nx * ny * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<nx * ny * 3; i++) indii[i]=-1 ;

	for (i=0; i<ny; i++) {
		latval = startlat - i * gspace ;
		for (j=0; j<nx ; j++) {
			lonval = startlon + j * gspace ;
			gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

			if (firstgrid) {
				listgrids[ngrids++]=gridnum ;
				firstgrid=0 ;
				indii [i * nx * 3 + j * 3] = 0 ;
			}
			else {
				foundgrid=0 ;
				// not a new grid
				for (igrid=0; igrid<ngrids; igrid++) {
					if (gridnum == listgrids[igrid]) {
						indii [i * nx * 3 + j * 3] = igrid ;
						foundgrid=1 ;
						break ;
					}
				}
				// if a new grid - unique number
				if (foundgrid==0) {
					indii [i * nx * 3 + j * 3] = ngrids ;
					listgrids[ngrids++] = gridnum ;
				}
			}
					
			//indii[i*nx*3+j*3]=gridnum ;
			indii[i*nx*3+j*3+1]=xloc ;
			indii[i*nx*3+j*3+2]=yloc ;
		}
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L * ngrids] ;
	float *indata_stdv = new float [1200L * 1200L * ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		// note that these are bsq files, mean, then stdev
		//sprintf (infile, "/local/worldbase/021km/orig_rad221_low_%03d_03_%02d_%02d.bsq\0",
                        sprintf (infile, "/local/worldbase/021km/orig/rad_%03d_03_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_03_%02d_%02d.bsq\0",
			mday, irow, icol) ;
		fin = fopen (infile, "r") ;
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdv[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < nx * ny ; isamp++) {
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		basegrid[isamp*2] = indata [inloc] ;
		basegrid[isamp*2+1] = indata_stdv [inloc] ;
		//mnstdev[isamp*2+1] = indata_stdv [inloc] ;
	}
/*
	FILE *fout =fopen ("testfile", "w") ;
	fwrite (basegrid, 4,  2*nx*ny, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdv ;
	delete [] indii ;



}
	
/* given a lat and a lon array, this function returns the 4 micron (221) mean and stand deviation 
 * the modisflag determines aqua or terra, night or day, (A_night T_day A_day T_night)
 * lat lon in decimal degrees and -180 to 180E longitude
 */
int sinu_1km::get_mn_stdev_array_b32 (float *lat, float *lon, int npix, int modisflag, float *mnstddev) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [npix * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<npix* 3; i++) indii[i]=-1 ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		latval = lat[i] ;
		lonval = lon[i] ;
		gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

		if (firstgrid) {
			listgrids[ngrids++]=gridnum ;
			firstgrid=0 ;
			indii [i*3] = 0 ;
		}
		else {
			foundgrid=0 ;
				// not a new grid
			for (igrid=0; igrid<ngrids; igrid++) {
				if (gridnum == listgrids[igrid]) {
					indii [i*3] = igrid ;
					foundgrid=1 ;
					break ;
				}
			}
				// if a new grid - unique number
			if (foundgrid==0) {
				indii [i*3] = ngrids ;
				listgrids[ngrids++] = gridnum ;
			}
		}
					
			//indii[i*nx*3+j*3]=gridnum ;
		indii[i*3+1]=xloc ;
		indii[i*3+2]=yloc ;
	
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L *  ngrids] ;
	float *indata_stdev = new float [1200L * 1200L *  ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		cout << listgrids[i] << "  : " << irow << "  :  " << icol << endl ;
		sprintf (infile, "/local/worldbase/021km/orig_rad32_low_%03d_0%1d_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_0%1d_%02d_%02d.bsq\0",
			mday, modisflag, irow, icol) ;
		fin = fopen (infile, "r") ;
		if (fin ==NULL) {
			cout << "UHOH could not read " << endl ;
			cout << infile << endl ;
			return(-1) ;
		}
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdev[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < npix ; isamp++) {
		if (badpix[isamp]) {
			mnstddev[isamp*2] = 0 ;
			mnstddev[isamp*2+1] = 10 ;
			continue ;
		}
			
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		mnstddev[isamp*2] = indata [inloc] ;
		mnstddev[isamp*2+1] = indata_stdev [inloc] ;
	}

/*
	FILE *fout =fopen ("testfile_mnsdv", "w") ;
	fwrite (mnstddev, 4,  npix * 2, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdev ;
	delete [] indii ;

	return (1) ;



}

/* given a lat and a lon array, this function returns the 4 micron (221) mean and stand deviation 
 * the modisflag determines aqua or terra, night or day, (A_night T_day A_day T_night)
 * lat lon in decimal degrees and -180 to 180E longitude
 */
int sinu_1km::get_mn_stdev_array_nti(float *lat, float *lon, int npix, int modisflag, float *mnstddev) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [npix * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<npix* 3; i++) indii[i]=-1 ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		latval = lat[i] ;
		lonval = lon[i] ;
		gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

		if (firstgrid) {
			listgrids[ngrids++]=gridnum ;
			firstgrid=0 ;
			indii [i*3] = 0 ;
		}
		else {
			foundgrid=0 ;
				// not a new grid
			for (igrid=0; igrid<ngrids; igrid++) {
				if (gridnum == listgrids[igrid]) {
					indii [i*3] = igrid ;
					foundgrid=1 ;
					break ;
				}
			}
				// if a new grid - unique number
			if (foundgrid==0) {
				indii [i*3] = ngrids ;
				listgrids[ngrids++] = gridnum ;
			}
		}
					
			//indii[i*nx*3+j*3]=gridnum ;
		indii[i*3+1]=xloc ;
		indii[i*3+2]=yloc ;
	
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L *  ngrids] ;
	float *indata_stdev = new float [1200L * 1200L *  ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		cout << listgrids[i] << "  : " << irow << "  :  " << icol << endl ;
		sprintf (infile, "/local/worldbase/021km/orig_radnti_low_%03d_0%1d_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_0%1d_%02d_%02d.bsq\0",
			mday, modisflag, irow, icol) ;
		fin = fopen (infile, "r") ;
		if (fin ==NULL) {
			cout << "UHOH could not read " << endl ;
			cout << infile << endl ;
			return(-1) ;
		}
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdev[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < npix ; isamp++) {
		if (badpix[isamp]) {
			mnstddev[isamp*2] = 0 ;
			mnstddev[isamp*2+1] = 10 ;
			continue ;
		}
			
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		mnstddev[isamp*2] = indata [inloc] ;
		mnstddev[isamp*2+1] = indata_stdev [inloc] ;
	}

/*
	FILE *fout =fopen ("testfile_mnsdv", "w") ;
	fwrite (mnstddev, 4,  npix * 2, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdev ;
	delete [] indii ;

	return (1) ;



}
/* given a lat and a lon array, this function returns the 4 micron (221) mean and stand deviation 
 * the modisflag determines aqua or terra, night or day, (A_night T_day A_day T_night)
 * lat lon in decimal degrees and -180 to 180E longitude
 */
int sinu_1km::get_mn_stdev_array_old (float *lat, float *lon, int npix, int modisflag, float *mnstddev) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [npix * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<npix* 3; i++) indii[i]=-1 ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		latval = lat[i] ;
		lonval = lon[i] ;
		gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

		if (firstgrid) {
			listgrids[ngrids++]=gridnum ;
			firstgrid=0 ;
			indii [i*3] = 0 ;
		}
		else {
			foundgrid=0 ;
				// not a new grid
			for (igrid=0; igrid<ngrids; igrid++) {
				if (gridnum == listgrids[igrid]) {
					indii [i*3] = igrid ;
					foundgrid=1 ;
					break ;
				}
			}
				// if a new grid - unique number
			if (foundgrid==0) {
				indii [i*3] = ngrids ;
				listgrids[ngrids++] = gridnum ;
			}
		}
					
			//indii[i*nx*3+j*3]=gridnum ;
		indii[i*3+1]=xloc ;
		indii[i*3+2]=yloc ;
	
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L *  ngrids] ;
	float *indata_stdev = new float [1200L * 1200L *  ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		cout << listgrids[i] << "  : " << irow << "  :  " << icol << endl ;
		sprintf (infile, "/local/worldbase/021km/orig_rad221_low_%03d_0%1d_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_0%1d_%02d_%02d.bsq\0",
			mday, modisflag, irow, icol) ;
		fin = fopen (infile, "r") ;
		if (fin ==NULL) {
			cout << "UHOH could not read " << endl ;
			cout << infile << endl ;
			return(-1) ;
		}
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdev[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < npix ; isamp++) {
		if (badpix[isamp]) {
			mnstddev[isamp*2] = 0 ;
			mnstddev[isamp*2+1] = 10 ;
			continue ;
		}
			
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		mnstddev[isamp*2] = indata [inloc] ;
		mnstddev[isamp*2+1] = indata_stdev [inloc] ;
	}

/*
	FILE *fout =fopen ("testfile_mnsdv", "w") ;
	fwrite (mnstddev, 4,  npix * 2, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdev ;
	delete [] indii ;


	return (1) ;

}

/* given a lat and a lon array, this function returns the 4 micron (221) mean and stand deviation 
 * the modisflag determines aqua or terra, night or day, (A_night T_day A_day T_night)
 * lat lon in decimal degrees and -180 to 180E longitude
 */
int sinu_1km::get_mn_stdev_array (float *lat, float *lon, int npix, int modisflag, float *mnstddev) {

	char infile[420] ;
	int thisngrid, i, j, igrid, ngrids=0, gridnum, imnum ;
	int listgrids [100], inloc, xloc, yloc, firstgrid=1, foundgrid ;
	int *indii = new int [npix * 3] ;
	float latval, lonval, xval, yval ;

	// initialize
	for (i=0; i<npix* 3; i++) indii[i]=-1 ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		latval = lat[i] ;
		lonval = lon[i] ;
		gridnum = get_gridnum (latval, lonval, &xloc, &yloc) ;  

		if (firstgrid) {
			listgrids[ngrids++]=gridnum ;
			firstgrid=0 ;
			indii [i*3] = 0 ;
		}
		else {
			foundgrid=0 ;
				// not a new grid
			for (igrid=0; igrid<ngrids; igrid++) {
				if (gridnum == listgrids[igrid]) {
					indii [i*3] = igrid ;
					foundgrid=1 ;
					break ;
				}
			}
				// if a new grid - unique number
			if (foundgrid==0) {
				indii [i*3] = ngrids ;
				listgrids[ngrids++] = gridnum ;
			}
		}
					
			//indii[i*nx*3+j*3]=gridnum ;
		indii[i*3+1]=xloc ;
		indii[i*3+2]=yloc ;
	
	}

	cout << "number of grids is " << ngrids << endl ;
	for (i=0; i< ngrids; i++) {
		cout << i << "  :  " << listgrids[i] << endl ;
		
	}

	float *indata = new float [1200L * 1200L *  ngrids] ;
	float *indata_stdev = new float [1200L * 1200L *  ngrids] ;
	FILE *fin ;
	int irow, icol , isamp;
	for (i=0; i<ngrids; i++) {
		irow = listgrids[i] / 36;
		icol = listgrids[i] - irow * 36 ;
		cout << listgrids[i] << "  : " << irow << "  :  " << icol << endl ;
		sprintf (infile, "/local/worldbase/021km/orig/rad_low_%03d_0%1d_%02d_%02d.bsq\0",
		//sprintf (infile, "/local/worldbase/021km/new_rad221_low_%03d_0%1d_%02d_%02d.bsq\0",
			mday, modisflag, irow, icol) ;
		fin = fopen (infile, "r") ;
		if (fin ==NULL) {
			cout << "UHOH could not read " << endl ;
			cout << infile << endl ;
			return(-1) ;
		}
		fread (&indata[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fread (&indata_stdev[i * 1200L * 1200L], 4, 1200L * 1200L, fin);
		fclose (fin) ;
	}
	
	for (isamp =0 ; isamp < npix ; isamp++) {
		if (badpix[isamp]) {
			mnstddev[isamp*2] = 0 ;
			mnstddev[isamp*2+1] = 10 ;
			continue ;
		}
			
		gridnum = indii[isamp*3] ;
		if (gridnum == 1) {
			int hg = 0 ;
			hg*=2 ;
		}
		xloc = indii[isamp*3+1] ;
		yloc = indii[isamp*3+2] ;
		inloc = gridnum * 1200L * 1200L + yloc * 1200L + xloc ;
		mnstddev[isamp*2] = indata [inloc] ;
		mnstddev[isamp*2+1] = indata_stdev [inloc] ;
	}

/*
	FILE *fout =fopen ("testfile_mnsdv", "w") ;
	fwrite (mnstddev, 4,  npix * 2, fout) ;
	//fwrite (indii, 4, 3 * nx *  ny, fout) ;
	fclose(fout) ;
	*/
		
	delete [] indata ;	
	delete [] indata_stdev ;
	delete [] indii ;


	return (1) ;

}


/**
 * allows the user to set the radius and central meridian if different than the default.
 * @param radius spherical radius of the planet - default is 6371007.181 m
 * @param centmeridian central meridian specifying the x zero location, long to 
 * the west is -ve, to the east is +ve -default is zero
 */
void sinu_1km::set_parameters (float radius, float centmeridian) {
    Radius = radius ;
    cent_merid = centmeridian ;
}

// latlon in degrees, gets converted to radians in function 
void sinu_1km::latlon_to_xy(float lat, float lon, float* x, float* y) {
    double latd = lat ;
    double lond = lon ;
    *x = Radius * (lond-cent_merid) * mDegToRad * cos (latd * mDegToRad) ;
    *y = Radius * latd * mDegToRad ;
}

void sinu_1km::xy_to_latlon (float x, float y, float *lat, float *lon) {
    *lat = y / Radius /mDegToRad ;
    *lon = cent_merid +( x / (Radius * cos(*lat * mDegToRad)))/mDegToRad ;
}

sinu_1km::sinu_1km(const sinu_1km& orig) {
}

sinu_1km::~sinu_1km() {
}
 
/**
 * 
 * @param mon       int - month to fill the grid with, should add something for which satellite, etc
 * @param aq_terra_flag int - aqua (0-night 2-day) terra (3-night 1-day)
 * @param ulc_lat   float - northernmost latitude
 * @param ulc_lon   float - westernmost longitude
 * @param gridspace float - grid spacing in decimal degrees 
 * @param nx        int - number of samples in output grid
 * @param ny        int - number of lines in output grid
 * @param outarr    float * - 8 banded nx * ny float array, memory allocated before calling function
 */
void sinu_1km::fillGrid (int mon, int aq_terra_flag, float ulc_lat, float ulc_lon, float gridspace, int nx, int ny, float *outarr) {
    
    int i, j, ixloc, iyloc, itype, npix_grid, iretn ;
    float lat, lon, xloc, yloc ;
    char bfile [420] ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    
    npix_grid = nx * ny ;
    
    // the ifil will depend upon the band, and the month to be extracted from the ave,stdev sinusoidal files.
    char const *outstr[] = {"rad221", "rad32", "nti"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;

    

    for (itype = 0; itype < 3; itype++) {
        sprintf(basef, "%s_%03d_%02d.bsq", outstr[itype], mon_days[mon], aq_terra_flag) ;
        sprintf(bfile, "/local/worldbase/02ssh/v2/%s_%03d_%02d.bsq\0", outstr[itype], mon_days[mon], aq_terra_flag);
        iretn = m02ssh_open (outstr[itype],mon_days[mon], modis_period(aq_terra_flag)) ;
        if (iretn <0) {
            cout << "Problem opening ... " << bfile << endl ;
        }
        cout << "opening baseline file of : " <<bfile << endl ;
        //FILE *fil = fopen(bfile, "r");
        //fread(aldat, 8, 4320L * 8640, fil);
        //fclose(fil);
        lst_val.jday = mon_days[mon] ;
        for (i = 0; i < ny; i++) {
            lat = ulc_lat - gridspace *  (i + 0.5) ;
            mpix.clear() ;
            lst_val.lat = RADOF (lat) ;
            
            for (j = 0; j < nx; j++) {
                lon = ulc_lon + gridspace * (j + 0.5)  ;
                lst_val.lon = RADOF (lon) ;
                mpix.push_back (lst_val) ;
                //latlon_to_xy (lat, lon, &xloc, &yloc) ;
                //iyloc = int (yloc) ;
                //ixloc = int (xloc) ;
                //outarr [2 * itype* npix_grid + i * nx + j] = 
                //outarr [2 * itype * npix_grid + npix_grid + i * nx + j] = aldat [8640L * 4320 + iyloc * 8640L + ixloc];
            }
            m02ssh_read (mpix) ;
            for (j=0; j<nx; j++) {
                outarr [2 * itype * npix_grid + i * nx + j] = mpix[j].mean ;
                outarr [2 * itype * npix_grid + npix_grid + i * nx + j] = mpix[j].std ;
            }
        }
    }
    /*
    FILE *fout = fopen ("/home/harold/aldat", "w") ;
    fwrite (outarr, 8, nx * ny, fout) ;
    fclose (fout) ;
    */
    
}



    


/**
 * getNTIValues will load up the outarr with nti values for the /local/worldbase/02ssh directory
 * the month, and correct aqua/terra flag are arguments to specify which file to choose from 
 * @param mon
 * @param aq_terra_flag
 * @param lat
 * @param lon
 * @param npix
 * @param outarr
 */
/* TEMPORARILY returning band 22values in */
void sinu_1km::getNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *outnti, bool openFlag) {
    
    int i ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    mpix.clear() ;
    
    
    // the ifil will depend upon the band, and the month to be extracted from the ave,stdev sinusoidal files.
    char const *outstr[] = {"rad221", "rad32", "nti"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;

    
    m02ssh_open(outstr[2], mon_days[mon], modis_period(aq_terra_flag));
    cout << endl;
    cout << endl;
    cout << "Opening for : " << outstr[2] << endl;
    cout << "Days :" << mon_days[mon] << endl;
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
     
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        outnti[i*2]= mpix[i].mean ;
        outnti[i*2+1]= mpix[i].std ;
    }
    
    
    // also get the band 22 values
    m02ssh_open(outstr[0], mon_days[mon], modis_period(aq_terra_flag));
    
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
    /* 
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
     */
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        out22[i*2]= mpix[i].mean ;
        out22[i*2+1]= mpix[i].std ;
    }
    

}



/**
 * getNTIValues will load up the outarr with nti values for the /local/worldbase/02ssh directory
 * the month, and correct aqua/terra flag are arguments to specify which file to choose from 
 * @param mon
 * @param aq_terra_flag
 * @param lat
 * @param lon
 * @param npix
 * @param outarr
 * @param outmax from the count file, band 21 max is band 2 of the count file
 */
void sinu_1km::getAllValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outnti, bool openFlag) {
    
    int i ;
    char countpath[240] ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    mpix.clear() ;
    
    
    // the ifil will depend upon the band, and the month to be extracted from the ave,stdev sinusoidal files.
    char const *outstr[] = {"rad221", "rad32", "nti", "count"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;

    
    m02ssh_open(outstr[2], mon_days[mon], modis_period(aq_terra_flag));
    cout << endl;
    cout << endl;
    cout << "Opening for : " << outstr[2] << endl;
    cout << "Days :" << mon_days[mon] << endl;
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
     
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        outnti[i*2]= mpix[i].mean ;
        outnti[i*2+1]= mpix[i].std ;
    }
    
    
    // also get the band 32 values
    m02ssh_open(outstr[1], mon_days[mon], modis_period(aq_terra_flag));
    
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
    
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
    
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        out32[i*2]= mpix[i].mean ;
        out32[i*2+1]= mpix[i].std ;
    }
    
    // also get the band 32 values
    m02ssh_open(outstr[0], mon_days[mon], modis_period(aq_terra_flag));
    
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
    /* 
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
     */
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        out22[i*2]= mpix[i].mean ;
        out22[i*2+1]= mpix[i].std ;
    }
    
    //sprintf (countstring, "/local/worldbase/02ssh/count_%03d_%02d")
//    sprintf(countpath, "/local/worldbase/02ssh/count_%03u_%02hhu.bsq", mon_days[mon], static_cast<uint8_t>(aq_terra_flag));
//    FILE *fin = fopen (countpath, "r") ;
//    unsigned short *cdata = new unsigned short [4320L * 8640L] ;
//    int x, y ;
//    fseek (fin, 2 * 4320L * 8640, SEEK_SET) ;
//    fread (cdata, 2, 4320L * 8640, fin) ;
//    for (i=0; i< mpix.size(); i++){
//        y = mpix[i].vindex ;
//        x = mpix[i].hindex ;
//        outmax21[i] = cdata[y * 8640 + x] / 1000. ;
//    }
//    fclose (fin) ;
//    delete [] cdata ;
    
}

void sinu_1km::getMaxNTIValues (int mon, int aq_terra_flag, float *lat, float *lon, int npix, float *out22, float *out32, float *outmaxNTI, bool openFlag) {
    
    int i ;
    char countpath[240] ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    mpix.clear() ;
    
    
    // the ifil will depend upon the band, and the month to be extracted from the ave,stdev sinusoidal files.
    char const *outstr[] = {"rad221", "rad32", "nti", "count"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;

/*    
    m02ssh_open(outstr[3], mon_days[mon], modis_period(aq_terra_flag));
    cout << endl;
    cout << endl;
    cout << "Opening for : " << outstr[2] << endl;
    cout << "Days :" << mon_days[mon] << endl;
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
     
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        outmax21[i*2]= mpix[i].mean ;
        outmax21[i*2+1]= mpix[i].std ;
    }
  */  
    
    // also get the band 32 values 
    m02ssh_open(outstr[1], mon_days[mon], modis_period(aq_terra_flag));
    
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
    
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
    
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        out32[i*2]= mpix[i].mean ;
        out32[i*2+1]= mpix[i].std ;
    }
    
    // also get the band 32 values
    m02ssh_open(outstr[0], mon_days[mon], modis_period(aq_terra_flag));
    
        //cout << "Modis type :  " << modis_period(aq_terra_flag) << endl ;
    /* 
    for (i=0; i<npix; i++) {
        lst_val.lat = RADOF (lat[i]) ;
        lst_val.lon = RADOF (lon[i]) ;
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period(aq_terra_flag) ;
        mpix.push_back (lst_val) ;
    }
     */
    m02ssh_read (mpix) ;
    
    for (i=0; i<npix; i++) {
        out22[i*2]= mpix[i].mean ;
        out22[i*2+1]= mpix[i].std ;
    }
    
    //sprintf (countstring, "/local/worldbase/02ssh/count_%03d_%02d")
    sprintf(countpath, "/local/worldbase/02ssh/count_%03u_%02hhu.bsq", mon_days[mon], static_cast<uint8_t>(aq_terra_flag));
    FILE *fin = fopen (countpath, "r") ;
    unsigned short *cdata = new unsigned short [4320L * 8640L] ;
    int x, y ;
    float val ;
    fseek (fin, 2 * 4320L * 8640, SEEK_SET) ;
    fread (cdata, 2, 4320L * 8640, fin) ;
    for (i=0; i< mpix.size(); i++){
        y = mpix[i].vindex ;
        x = mpix[i].hindex ;
        val = cdata[y * 8640 + x] / 1000. ;
        outmaxNTI[i] = (val - out32[i*2]) / (val + out32[i*2]) ;
    }
    fclose (fin) ;
    delete [] cdata ;
    
}

void sinu_1km::getCorrespondingValues(vector<int> alinds, int mon, int aq_terra_flag,  float *lat, float *lon, float *vals) {
    char bfile [240] ;
    int i, ind, itype, num_alerts, iretn ;
    char const *outstr[] = {"rad221" , "rad32", "nti"} ;
    int mon_days[] = {1, 32,62, 92,123, 153,184,214,245,275,305,336} ;
    int ixloc, iyloc ;
    float xloc, yloc ;
    m02ssh_coordvalue lst_val ;
    vector <m02ssh_coordvalue> mpix ;
    
    num_alerts = alinds.size() ;
    
    if (num_alerts <= 0) return ;
    for (itype = 0; itype < 3; itype++) {
        mpix.clear() ;
        sprintf(basef, "%s_%03d_%02d.bsq", outstr[itype], mon_days[mon], aq_terra_flag) ;
        sprintf(bfile, "/local/worldbase/02ssh/v2/%s_%03d_%02d.bsq\0", outstr[itype], mon_days[mon], aq_terra_flag);
        iretn = m02ssh_open (outstr[itype],mon_days[mon], modis_period(aq_terra_flag)) ;
        if (iretn < 0)
            cout << "problem opening baseline file of : " <<bfile << endl ;
        
        lst_val.jday = mon_days[mon] ;
        lst_val.period = modis_period (aq_terra_flag) ;
        for (i = 0; i < num_alerts; i++) {
            ind = alinds.at(i);
            this->latlon_to_xy(lat[ind], lon[ind], &xloc, &yloc);
            lst_val.lat = RADOF (lat[ind]) ;
            lst_val.lon = RADOF (lon[ind]) ;
            mpix.push_back (lst_val) ;
            
            

        }
        m02ssh_read (mpix) ;
        for (i = 0; i < num_alerts; i++) {
            vals [i*6+2 * itype] = mpix[i].mean ;
            vals [i*6+2 * itype+1] = mpix[i].std ;
            
        }
        
    }
    
    
}
