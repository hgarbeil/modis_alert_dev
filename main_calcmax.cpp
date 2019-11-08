///
// This program - modis_calcmax - works off of the /local/worldbase/021km/orig rad221_highest and rad32_highest  files 
// This has the top 15 values (excluding >-.8nti alerts) encountered over a 10 year period for a pixel.
// Look at the series of values and find the max before a slope break up to an anomalous high
#include "sinu_1km.h"
#include <stdio.h>
#include <string.h>
#include <cstdint>
#include <iostream>
//#include "modis_hdf.h"
//#include "modis_process.h"
#include "alert.h"
#include "testime.h"
#include "median.h"

int main (int argc, char *argv[]) {
    
    char *tmp, *tmp1, flist[420], infile [420], intmp[120], outfile[420], prefix[420], rad32file[120] ;
    uint16_t *indat_22, *indat_32 ;
    float *histo_22, *histo_32, *diff, *temp_diff, *maxarr, val32, val22, maxval, maxTdiff ;
    float temp22, temp32 ;
    int i, j, ns, nl, npts, maxsub ;
	float inarr[9] ;
    FILE *fin, *fout ;	
                    if (argc < 3) {
		cout << "Usage: modis_calcmax filelist  outdir " << endl ;
		exit(-1) ;
	}

                    strcpy (flist, *++argv) ;
                    strcpy (prefix, *++argv) ;

                    alert *al = new alert() ;
                    int this_mday ;
                    ns = 1200 ;
                    nl = 1200 ;
                    npts = 15 ;

                    indat_22 = new uint16_t [ns * nl * npts] ; 
                    indat_32 = new uint16_t [ns * nl * npts] ; 
                    maxarr = new float [ns * nl*3] ;
                    histo_22 = new float [npts] ;
                    histo_32 = new float [npts] ;
                    diff = new float [npts] ;
                    temp_diff = new float [npts] ;
					float *tmparr=new float [nl * ns] ;

	

          FILE *flistfil = fopen (flist, "r") ;
          if (flistfil == NULL) {
                    cout << "could not open  " << flist << endl ;
          }
          while (!feof(flistfil)) {
                    fscanf (flistfil, "%s", intmp) ;
          
					strcpy (infile, "/local/worldbase/021km/orig/") ;
                    //strcpy (intmp, *++argv) ;
                    strcat (infile, intmp) ;
                    tmp = strstr (intmp, "t_") ;
                    strcpy (rad32file, "rad32_highes") ;
                    strcat (rad32file, tmp) ;
                    cout << "Rad 32 file is this... " << rad32file << endl ;

                    fin = fopen (infile, "r") ;
                    // make sure that input file opens
                    if (fin ==NULL) {
                        cout << "Could not open : " << infile << endl ;
                        fclose (fin) ;
                        return (-1) ;
                    }
                    fread (indat_22, 2, ns * nl * npts, fin) ;
                    fclose (fin) ;

                    strcpy (infile, "/local/worldbase/021km/orig/") ;
                    strcat (infile, rad32file) ;

                    fin = fopen (infile, "r") ;
                    if (fin ==NULL) {
                        cout << "Could not open : " << infile << endl ;
                        fclose (fin) ;
                        return (-1) ;
                    }
                    fread (indat_32, 2, ns * nl * npts, fin) ;
                    fclose (fin) ;


	//strcpy (prefix, *++argv) ;
                    strcpy (outfile, prefix) ;
                    strcat (outfile, "/") ;
                    strcat (outfile, intmp) ;
                    tmp = strstr (outfile, ".bsq") ;
                    strcpy (tmp, "_max") ;

                    // go through line by line
                    for (i=0; i<ns * nl; i++) {
                        // load up histogram
                        maxsub = 0 ;
                        maxval = 0 ;
                        maxTdiff = 0. ;
                        for (j=0; j<npts; j++) {
                            histo_22[j] = (float)indat_22[j * ns * nl+i] ;
                            histo_32[j] = (float)indat_32[j*ns*nl+i] ;
                            temp32= al->bb_radtotemp(12., histo_32[j]/1000.) ;
                            temp22= al->bb_radtotemp(4., histo_22[j]/1000.) ;
                            val22 = (float)indat_22[j * ns * nl+i] ;
                            val32 = (float)indat_32[j * ns * nl+i] ;
                            //temp_diff[j] = (val22 - val32)/(val22 + val32) ;
                            temp_diff[j] = temp22 - temp32 ;
                            // check if hotspot or a cloud, either case ignore....
                            if (temp_diff[j] > 20.) continue ;
                            
                             if (histo_22[j] > maxval) {
                                    maxval = histo_22[j] ;
                                    maxsub = j ;
                                    maxTdiff = temp_diff[j] ;
	                
                                 }
                            
                        }
                        // if never found a good pixel
                        if (maxsub==0) maxsub=0 ;
                        maxarr[i] = histo_22[maxsub]/1000. ;
                        maxarr[i + ns * nl] = histo_32[maxsub]/1000. ;
                        maxarr[i + ns * nl * 2] = maxsub ;
                        }

					// median filt
					/*****
					int ii, jj, isub, jsub, count ;
					for (i=0; i<ns * nl; i++) tmparr[i] = maxarr[i] ;

					for (i=0; i<nl; i++) {
						for (j=0; j<ns; j++) {
							count = 0 ;
							for (isub=-1; isub<=1; isub++) {
								ii = i + isub ;
								if (ii < 0) continue ;
								if (ii >= nl) continue ;
								for (jsub=-1; jsub<=1; jsub++) {
									jj = j + jsub ;
									if (jj < 0) continue ;
									if (jj >= ns) continue ;
									inarr[count] = tmparr[ii * ns + jj] ;
									count++  ;
								}
							}
							// now filter
							maxarr[i * ns + j] = fmedian (inarr, count) ;
						}
					}
					*****/
							



									



                    

                fout = fopen (outfile, "w") ;
				fwrite (maxarr, 4, 3*ns*nl, fout) ;
				fclose (fout) ;
          } // next file
					delete [] tmparr ;
                    delete [] indat_22 ;
                    delete [] indat_32 ;
                    delete [] histo_22 ;
                    delete [] histo_32 ;
                    delete [] maxarr ;
                    delete [] diff ;
                    delete [] temp_diff ;
}
