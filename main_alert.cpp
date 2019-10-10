/**
 * main_alert.cpp
 * Purpose : get the alert values for all pixels in the image
 * can put them out in a huge text file or in a resampled image file
 * @author Harold Garbeil HIGP/SOEST/UHM
 * @version 1.0 6/30/2018
 * 
 **/

 

#include <cstdlib>
#include "modis_hdf.h"
#include "modis_process.h"
#include "surftemp.h"
#include <boost/filesystem/path.hpp>

 
using namespace std;

/**
 * Main function of the modis_alert program.
 *  
 **/ 
int main(int argc, char** argv) {


      int count ;
     bool openFlag ;
    char outstr[800],  histofile[420], tfile [420], mfile[420], ofile[420], ofile1[420],flist[420], outpref[420], outputfile[420],  mstdv_file[420], modid[420] ;
    char alertfile[420], localertfile[420] ;
    char proclogfile [420] ;
    char *sufind, *sufind1 ;
    int minind ;
    float testtemp;
    float rlat, rlon, pvals[6];
    rlat = 19.333;
    rlon = -155.333;
    count = 0 ;
    int stuff [4], oldmonth = -1;
    modis_hdf *geom, *therm;
    modis_process *mproc;

    openFlag = true ;
    
    
    // get ASCII file with emiss, geom hdf files
    //strcpy (flist, *++argv) ;
    strcpy (flist, "/home/harold/workdir/alert_results/hawaii_aqua.txt") ;
    cout << "File list is " << flist << endl ;
    // get output path 
    strcpy (outpref, "/home/harold/workdir/alert_results/modis_dev/hawaii/aqua_night") ;
    strcpy (proclogfile, outpref) ;
    strcat (proclogfile, "/proclog.txt\0") ;
    strcpy (histofile,  "/home/harold/workdir/alert_results/modis_dev/hawaii/aqua_night/aqua_night_histo") ;
    strcpy (alertfile, "/home/harold/workdir/alert_results/modis_dev/hawaii/aqua_night/aqua_night_alert.txt") ;
    strcpy (localertfile, "/home/harold/workdir/alert_results/modis_dev/hawaii/aqua_night/aqua_night_localert.txt") ;
    cout << "outpref is " << outpref << endl ;
    // and the log file
    
   
    
    
    cout << "File list is : " << flist << endl ;
    
    FILE *fin = fopen(flist, "r");
   
    if (fin == NULL) {
        cout << " Could not open file list :  " << flist << endl;
        exit(-1);
    }
    // controlling process class
    mproc = new modis_process();
    mproc->set_proclog (proclogfile) ;
    
    
    // go through each pair of MODIS MOD021KM and corresponding MOD03 files
    // calculating NTI values for each pixel in the 1354 x 2030 original array
    // as well as resampling radiance bands and NTI values to the predetermined
    // Hawaii ROI
    while (!feof(fin) && count < 12000) {
        fscanf(fin, "%s %s", tfile, mfile);
        
        sufind = strstr (mfile, "A2014") ;
        strncpy (modid, &mfile[sufind-mfile], 13) ;
        modid[13] = '\0' ;
        therm = new modis_hdf(tfile);
        therm->get_file_name (tfile, ofile) ;
        strcpy (ofile1, ofile) ;
        
        sufind = strstr (ofile, ".hdf") ;
        
        
        strncpy (sufind, "_out", 4) ;
        strcpy (outputfile, outpref) ;
        strcat (outputfile, "/") ;
        
        strcat (outputfile, ofile) ;
        
        cout << "Output cube file will be " << outputfile << endl ;
        
        // alert text file
        sufind = strstr (ofile1, ".hdf") ;
        strncpy (sufind, "_alert.txt\0", 11) ;
//        strcpy (alertfile, outpref) ;
//        strcat (alertfile, "/") ;
//        strcat (alertfile, ofile1) ;
//        cout << "Alert file will be " << alertfile << endl ;
        sufind = strstr (ofile1, "_alert.txt") ;
        strncpy (sufind, "_localert.txt\0", 14) ;
//        strcpy (localertfile, outpref) ;
//        strcat (localertfile,"/") ;
//        strcat (localertfile, ofile1) ;
        
        
        
         
        
        therm->get_date_period(tfile, stuff);
        cout << "HDF file month is : " << stuff [1] << endl ;
        geom = new modis_hdf(mfile);
        if (!geom->geom_status || geom->dayflag || therm->readstatus<0) {
            delete geom ;
            delete therm ;
            continue ; 
        }
        mproc->set_month(stuff[1]-1);
        
        //if (mproc->mon < 10) continue ;
        
        mproc->set_modis_hdfs(geom, therm);
        mproc->set_prefix (modid) ;
        if (count==0)
            //mproc->set_bounds (14.7563,-91.5534, .008, 400,400) ;
            //mproc->set_bounds (13.61, 40.67, .008 ,700,700) ;
            // hawaii
            mproc->set_bounds (20.3, -156, .008,1200, 1200) ;
            //mproc->set_bounds (-23.33, -67.75, .008, 400, 400) ;
            // iceland
            //mproc->set_bounds (64.63, -17.53,.008,700, 700) ;
            //nyiragongo
            //mproc->set_bounds (-1.55, 29.25, .008, 700, 700) ;
        
        
        if (mproc->mon !=oldmonth ) {
            mproc->extract_from_baseline_file() ;
            oldmonth = mproc->mon ;
            openFlag = true ;
        }
        cout << "*******************************" << endl ;
        cout << "Processing image number : " << count << endl ;
        
        // calc alert is the full scene alert calculation- this fills the alinds vector
        // a second alert algorithm is added to the calculation 
        // in this, we get the global avg nti and stdv and then compare it to the scene nti value, if
        // the scene value exceeds, then we mark it as an alert....
        mproc->calc_alert_bb(alertfile, localertfile) ;
        //mproc->calc_alert(alertfile, localertfile, openFlag) ;
        openFlag = false ;
        cout<< "calc_alert done " << endl ;
        
        // fill the grid for output image file, only for inspection purpose, calc_alert goes through entire image
        //cout << "grid process " << endl ;
        if (count < 50)
        mproc->process() ;
        
        //mproc->procalert() ;
        
        //mproc->get_nearest_pixel(rlat, rlon, &minind, pvals);
        cout << "Dayflag is " << geom->dayflag << endl;
        
        //strcpy (mstdv_file, "/local/worldbase/02ssh/rad32_001_00.bsq") ;
        // get the alert from MODIS file
        
        
        //mproc->zscore() ;
        mproc->write_output (outputfile) ;
        cout << flush;
        
        delete geom;
        delete therm; 
        
        count++ ;

    }
    

    fclose(fin);
    FILE *fhisto = fopen (histofile, "w") ;
    fwrite (mproc->histoarr, 2, mproc->ny_grid * mproc->nx_grid * 6, fhisto) ;
    fclose (fhisto) ;

    delete mproc;
   // fclose(fout);
    return 0;
}

