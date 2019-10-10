/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: harold
 *
 * Created on May 1, 2018, 12:09 PM
 */

#include <cstdlib>
#include "modis_hdf.h"
#include "modis_process.h"
#include "surftemp.h"


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {


	vector<lst_coord> mc ;
	lst_coord tcoord ;
	char outstr[800], tfile [420], mfile[420], ofile[420], flist[420] ;
	int minind ;
	float testtemp ;
	float rlat, rlon,  pvals[6] ;
	rlat = 19.333 ;
	rlon = -155.333 ;
	
	int stuff [4] ;
	modis_hdf *geom, *therm ;
	modis_process *mproc ;
	strcpy (flist, argv[1]) ;
/*
	strcpy (tfile, *++argv) ;
	strcpy (mfile, *++argv) ;
	strcpy (ofile, *++argv) ;
	geom = new modis_hdf (mfile) ;
*/
	FILE *fout = fopen ("nightfiles.txt", "w") ;
	FILE *fin = fopen (flist, "r") ;
	if (fin == NULL) {
		cout << " Could not open " << flist << endl ;
		exit(-1) ;
	}
	mproc = new modis_process () ;
	while (!feof(fin)) {
		fscanf (fin, "%s %s",tfile, mfile) ;

	
		therm = new modis_hdf (tfile) ;
		therm->get_date_period(tfile, stuff) ; 
		geom = new modis_hdf (mfile) ;
		mproc->set_modis_hdfs (geom, therm) ;
		mproc->set_month (stuff[1]) ;
	//mproc->procalert() ;
	//mproc->write_output ("output.dat") ;
		mproc->get_nearest_pixel (rlat, rlon, &minind, pvals) ;
		cout << "Dayflag is " << geom->dayflag << endl ;
		if (!geom->dayflag && pvals[0]<1200.) 
		{
			sprintf (outstr, 
				"%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f", 
				tfile, stuff[0], stuff[1], stuff[2], stuff[3], 
				pvals[0], pvals[5], pvals[1], pvals[2], pvals[3],pvals[4]) ;	
			fprintf (fout, "%s\r\n", outstr) ;
			fflush (fout) ;
		}
		cout << flush ;
		delete geom ;
		delete therm ;

	}
    delete mproc ;

    fclose (fin) ;

    fclose (fout) ;
    return 0;
}

