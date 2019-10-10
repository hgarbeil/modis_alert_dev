#include "modis_process.h"
#include <math.h>
#include "surftemp.h"
#include "sinuProjection.h"

modis_process::modis_process() {

    distarr = 0L;
    bandsarr = 0L;
    temp = 0L;
    alerts = 0L;
    nl = 2030;
    ns = 1354;
    gridspace = 0.1;
    day_limit = -.6;
    night_limit = -.8;
    alerts = new float [ns * nl];
    newalerts = new float [ns * nl] ;
    sp = new sinuProjection () ;
    glob_mnstdev = 0L ;
    basevals_full = 0l ;
    basevals_full_new = 0l ;
    strcpy (proclog, "") ;
    alertFlag = false ;


}

modis_process::~modis_process() {

    delete [] bandsarr;
    delete [] distarr;
    delete [] alerts;
    delete [] newalerts ;
    delete [] histoarr ;
    if (temp) delete [] temp;
    if (basevals_full) delete[] basevals_full ;
    delete [] glob_mnstdev ;
    delete sp ;
}

void modis_process::set_modis_hdfs(modis_hdf *g, modis_hdf *th) {
    geom = g;
    therm = th;
    alertFlag = false ;
}

void modis_process::set_bounds (float cent_lat, float cent_lon, float gspace, int ny, int nx) {
    int i, npix_grid ;
    ny_grid = ny ;
    nx_grid = nx ;
    npix_grid = ny_grid * nx_grid ;
    gridspace = gspace ;
    startlat = cent_lat + ny/2 * gspace ;
    startlon = cent_lon - nx/2 * gspace ;
    endlat = startlat - (ny_grid-1) * gridspace ;
    endlon = startlon + (nx_grid -1) * gridspace ;
    if (distarr) {
        delete [] distarr;
        delete [] bandsarr;
        delete [] glob_mnstdev ;
    }
    distarr = new float [npix_grid];
    bandsarr = new float [7 * npix_grid];
    glob_mnstdev = new float [6 * npix_grid] ;
    histoarr = new unsigned short [6 * npix_grid] ;
    // alloc memory for the baseline data b221, b32, may do all b6 and ndti
    
    for (i = 0; i < npix_grid; i++) {
        distarr [i] = 9999.;
        histoarr[i+npix_grid] = 0 ;
        histoarr[i+npix_grid*2] = 0 ;
        histoarr[i+npix_grid*3] = 0 ;
        histoarr[i+npix_grid*4] = 0 ;
        histoarr[i+npix_grid*5] = 0 ;
    }
    for (i = 0; i < npix_grid * 7; i++) {
        bandsarr [i] = -1.;
    }

    
}


void modis_process::set_bounds(float ulc_lat, float ulc_lon, float lrc_lat,
        float lrc_lon, float gspace) {

    int i, npix_grid ;


    gridspace = gspace;

    ny_grid = int (((ulc_lat - lrc_lat) / gridspace) + 1);
    nx_grid = int (((lrc_lon - ulc_lon) / gridspace) + 1);
    npix_grid = nx_grid * ny_grid ;
    startlat = ulc_lat;
    startlon = ulc_lon;
    endlat = lrc_lat;
    endlon = lrc_lon;

    cout << "Number of lines is : " << ny_grid << endl;
    cout << "Number of samples is : " << nx_grid << endl;

    if (distarr) {
        delete [] distarr;
        delete [] bandsarr;
        delete [] glob_mnstdev ;
    }
    distarr = new float [npix_grid];
    bandsarr = new float [7 * npix_grid];
    histoarr = new unsigned short [5 * npix_grid] ;
    glob_mnstdev = new float [6 * npix_grid] ;
    // alloc memory for the baseline data b221, b32, may do all b6 and ndti
    
    for (i = 0; i < npix_grid; i++) {
        distarr [i] = 9999.;
        histoarr[i] = 0 ;
        histoarr[i+npix_grid] = 0 ;
        histoarr[i+npix_grid*2] = 0 ;
        histoarr[i+npix_grid*3] = 0 ;
        histoarr[i+npix_grid*4] = 0 ;
        histoarr[i+npix_grid*5] = 0 ;
        
    }
    for (i = 0; i < npix_grid * 7; i++) {
        bandsarr [i] = -1.;
        
    }
    

    
}


void modis_process::set_month(int mon) {
    this->mon = mon;
}

void modis_process::set_prefix (char *pref){
    strcpy (modid, pref) ;
    
}

void modis_process::get_nearest_pixel(float lat, float lon, int *index, float *pixvals) {
    int i, minind;
    float plat, plon;
    double xdist, ydist, dist, mperdegree, mindist;

    mindist = 1.E9;
    mperdegree = 110000.;

    // go through the data set lat and lon arrays :finding pixel closest to the ref pixel
    // defined by the function arguments (lat, lon)
    // the function will load up the index pointer with the closest pixel index and the
    // pixvals will be loaded with
    // [0] = mindist
    // [1] = b21 radiance
    // [2] = b22 b temp
    // [3] = b22 radiance
    // [4] = b32 radiance
    // [5] = solar zenith
    int npix = 1354L * 2030;
    for (i = 0; i < npix; i++) {
        ydist = (geom->latarr[i] - lat) * mperdegree;
        xdist = (geom->lonarr[i] - lon) * mperdegree;
        dist = sqrt(xdist * xdist + ydist * ydist);
        if (dist < mindist) {
            mindist = dist;
            minind = i;
        }
    }
    *index = minind;
    pixvals[0] = mindist;
    // b21
    pixvals[1] = therm->raddata_cal[minind];
    pixvals[2] = this->bb_radtotemp(3.992, pixvals[1]);
    // b22
    pixvals[3] = therm->raddata_cal[npix + minind];
    // b31
    pixvals[4] = therm->raddata_cal[2 * npix + minind];
    // solar zenith
    pixvals[5] = geom->solzen[minind] / 100.;
}
/*
void modis_process::procalert() {

    float plat, plon, dist_pixel;
    int i;
    size_t period;
    bool aflag, dflag;
    float b21val, b21temp;
    vector<lst_coord> mpix;
    lst_coord tcoord;
    int npix = 2030 * 1354;
    temp = new float [npix];

    // get the period 
    aflag = therm->aquaflag;
    dflag = therm->dayflag;


    // if aqua, daytime is 2, nighttime is 1
    if (aflag) {
        if (dflag)
            period = 2;
        else
            period = 0;
    }        // if terra, daytime is 1, nighttime is 3
    else {
        if (dflag)
            period = 1;
        else period = 3;
    }

    tcoord.period = lst_period(period);



    for (i = 0; i < npix; i++) {
        tcoord.lat = RADOF(geom->latarr[i]);
        tcoord.lon = RADOF(geom->lonarr[i]);
        plat = geom->latarr[i];
        plon = geom->lonarr[i];
        tcoord.month = mon;
        mpix.push_back(tcoord);

        b21val = therm->raddata_cal [i];
        // convert the b21 radiance to temperature
        b21temp = this->bb_radtotemp(3.992, b21val);
        temp[i] = b21temp;


    }


    // now get the mean and std dev
    lst_open();
    lst_read(mpix);
    float zval;
    for (i = 0; i < npix; i++) {
        tcoord = mpix.at(i);
        if (tcoord.mean < 200) {
            continue;
        }
        zval = (temp[i] - tcoord.mean) / tcoord.std;
        if (zval > 1) {
            //cout << "Mean temp is : " << tcoord.mean <<  "  "  << tcoord.std<< "  " << temp[i] << endl ;
        }
    }


    //delete [] temp ;	


}
*/
float modis_process::calcdist(float plat, float plon) {

    double xdist, ydist, tdist;
    float reflat, reflon;
    reflat = 19.33;
    reflon = -155.33;
    xdist = abs(reflat - plat) * 100000.;
    ydist = abs(reflon - plon) * 100000.;
    tdist = sqrt(xdist * xdist + ydist * ydist);
    return (float(tdist));
}

void modis_process::process() {

    int i, j, ib, snum, is, js, ixloc, iyloc, grid_x, grid_y, gridloc, npix_grid;
    float xloc, yloc, xdist, ydist, dist, latval, lonval;
    npix_grid = ny_grid * nx_grid;
    ;

    //if (!alertFlag) this->calc_alert() ;

    for (i = 0; i < npix_grid; i++) {
        distarr [i] = 9999.;
    }
    for (i=0; i< npix_grid * 4; i++) {
        bandsarr[i] = -1 ;
    }
    for (i = 0; i < 2030; i++) {
        for (j = 0; j < 1354; j++) {
             
            snum = i * 1354 + j;
            latval = geom->latarr[snum];
            lonval = geom->lonarr[snum];
            if (latval > startlat || latval < endlat)
                continue;
            if (lonval < startlon || lonval > endlon)
                continue;

           
            
            xloc = (lonval - startlon) / gridspace;
            yloc = (startlat - latval) / gridspace;
            ixloc = int(xloc + 0.5);
            iyloc = int(yloc + 0.5);

            for (is = -3; is <= 3; is++) {
            //is = 0 ;
                grid_y = iyloc + is;
                if (grid_y < 0) continue;
                if (grid_y >= ny_grid) continue;

                ydist = yloc - grid_y;
                for (js = -3; js <= 3; js++) {
                //js=0 ;
                    grid_x = js + ixloc;
                    if (grid_x < 0) continue;
                    if (grid_x >= nx_grid) continue;
                    gridloc = grid_y * nx_grid + grid_x;
                    xdist = xloc - grid_x;
                    dist = sqrt(xdist * xdist + ydist * ydist);
                    if (dist < distarr[gridloc]) {
                        distarr[gridloc] = dist;
                        //bandsarr[gridloc] = therm->refdata_cal[snum];
                        // resample the three bands of radiance to bandsarr
                        
                        for (ib = 0; ib < 3; ib++)
                            bandsarr[ib * npix_grid + gridloc] = therm->raddata_cal[ib * 2030L * 1354 + snum];
                        bandsarr[3*npix_grid+gridloc] = *(alerts+snum) ;

                        //bandsarr[4L * npix_grid + gridloc] = geom->solsens[snum];
                        //bandsarr[5L * npix_grid + gridloc] = geom->solsens[2030L * 1354 + snum];
                        //bandsarr[6L * npix_grid + gridloc] = geom->solsens[2 * 2030L * 1354 + snum];
                        //bandsarr[7L * npix_grid + gridloc] = geom->solsens[3 * 2030L * 1354 + snum];

                    }
                }
            }
        }
    }
}
/**
 * 
 * @param alertfile
 */



void modis_process::set_proclog (char *plog) {
    strcpy (proclog, plog) ;
}

/*
void modis_process::calc_alert_alice (char *alertfile, char *localertfile) {
    
    
}
*/

void modis_process::calc_alert_bb (char *alertfile, char *localertfile) {
    bool ingrid, stdflag ;
    int ixloc, iyloc, ind ;
    float zval, z_lon, z_lat, xloc, yloc;
    int i, j,  istep,  npix, npix_grid, num_alerts_new ;
    float std_alert, val31,btemp32, btemp21, btemp31, val21, val22, val32,  b12, b22threshold, planck22, rad_ratio, b22add ;
    int alert_count [6] ;
    float hotfrac []= {0.00015, .000175, .00020, .00025, .0005} ;
    char outstr[240] ;
    
    float *glob_21 = new float [nl * ns * 2] ;
    float *glob_32 = new float [nl * ns * 2] ;
    float *glob_nti = new float [nl * ns * 2] ;
    
    FILE *fproc = fopen (proclog, "a") ;
    for (i=0; i<6; i++)alert_count[i] = 0 ;
    char ntmp[20] ;
    int jday ;
    int utime ;
    
    strncpy (ntmp, &modid[5],3) ;
    ntmp[3]='\0';
    jday = atoi (ntmp) ;
    strncpy (ntmp, &modid[9],4) ;
    ntmp[4] ='\0' ;
    utime = atoi (ntmp) ;
    
    sp->getAllValues (mon, geom->aq_terra_flag, &geom->latarr[0], &geom->lonarr[0], ns*nl,glob_21, glob_32, glob_nti, true) ;
    b22add = this->bb_temptorad (3.95, 500+273.15) ;
    npix = ns * nl ;
    npix_grid = nx_grid * ny_grid ;
    
    // unique_alerts are those found with new method rather than std
    vector <int> new_alinds ;
    alinds.clear() ;
    
    
    
    for (int ij=0; ij<npix; ij++){
        //ij = i * this->ns + j ;
        val21 = therm->raddata_cal[ij];
        val22 = therm->raddata_cal[npix + ij];
        val31 = therm->raddata_cal[2 * npix + ij];
        val32 = therm->raddata_cal[3 * npix + ij];
        val32 = (val32 > glob_32[ij*2]) ? val32 : glob_32[ij*2] ;
        btemp32 = this->bb_radtotemp (12., val32) ;
        btemp31 = this->bb_radtotemp (11., val31) ;
        btemp32 = (btemp32 > btemp31) ? btemp32 : btemp31 ;
        
        planck22 = this->bb_temptorad (3.95, btemp32+10) ;
        
        z_lon = geom->lonarr[ij];
        z_lat = geom->latarr[ij];
        xloc = (z_lon - startlon) / gridspace;
        yloc = (startlat - z_lat) / gridspace;
        ixloc = int(xloc + 0.5);
        iyloc = int(yloc + 0.5);
        ingrid = true;
        if (ixloc < 0 || iyloc < 0) ingrid = false;
        if (iyloc >= ny_grid || ixloc >= nx_grid) ingrid = false;
        if (ingrid==true){
            int hg ;
            hg*=1 ;
        }
       
        
        // lpreprocess pixel radiance values     
        // check if val22 is saturated, then use val21
        // check if val32 is less than the mean value for 32
        //if (val32 < glob_32[ij*2]) continue ;
        //if either 21 | 22 or 32 < 0. continue 
        if (val22 > 2.0 || val22 < 0.001)
            val22 = val21;
        // also check to make sure 21 and 32 are good 
        if (val21 > 150. || val32 > 40.) {
            *(alerts + ij) = -2.;
            *(newalerts+ij) = -2. ;
            continue;
        }
        if (val21 < val22) val22 = val21 ;
        if (val22 <= 0. || val32 <= 0.) {
            *(alerts + ij) = -2.;
            *(newalerts+ij) = -2. ;
            continue;
        }
        //alval = (val22 - val32) / (val22 + val32);
        std_alert = (val22 - val32) / (val22 + val32);
        alerts[ij] = std_alert ;
        stdflag = false ;
        if (std_alert > -0.8) {
            alert_count[0] += 1 ;  
            stdflag = true ;
            if (ingrid) {
                    histoarr[iyloc*nx_grid+ixloc]+=1 ;
            }
        }
        for (istep=1; istep<6; istep++) {
            b22threshold = b22add * hotfrac[istep-1] + planck22 * 1 ;
            if (istep==2) {
                newalerts[ij] = val22 - b22threshold ;
                // write to alerts file if old alert is true but new is false
                if (stdflag && newalerts[ij]<0.){
                    alinds.push_back(ij) ;
                }
            }
            
            if (val22 > b22threshold) {
                alert_count[istep]++ ;
                // check if unique
                if (istep ==2 && stdflag == false){
                    new_alinds.push_back(ij) ;
                }
                if (ingrid) {
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid*istep]+=1 ;
                }
                if (istep == 2){
                    int hg = 0 ;
                    hg*= 2 ;
                }
            } else break  ;
        } 
    }
    
    sprintf (outstr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\r\n", jday, utime, alert_count[0], alert_count[1], alert_count[2], 
            alert_count[3], alert_count[4], alert_count[5]) ;
    fputs (outstr, fproc) ;
    
    FILE *foutnew = NULL ; 
    num_alerts_new = new_alinds.size() ;
    if (num_alerts_new > 0) {
    foutnew = fopen (localertfile, "a") ;
    for (i=0; i< num_alerts_new; i++) {
        
        ind = new_alinds.at(i) ;
        //zval = (alerts[ind] - basevals_full_new[6*i+4])/basevals_full_new[6*i+5] ;
        zval = (alerts[ind] - glob_nti[ind*2])/glob_nti[ind*2+1] ;
        //zval = alerts[ind] / glob_max21[ind] ;
         //if (zval < 1.1) continue ;
//        sprintf (outstr, "%d\t%s\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%3.1f\r\n", 
//                    ind, modid, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
//                    alerts[ind], zval) ;
        sprintf (outstr, "%d\t%d\t%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%7.3f\t%7.3f\t%7.3f\r\n", 
                    ind, jday, utime, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
                    alerts[ind], newalerts[ind], glob_21[ind*2], glob_32[ind*2]) ;
        //sprintf (outstr, "%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%6.3f\r\n", 
        //            ind, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
        //            alerts[ind], glob_nti[2*ind]) ;
        fputs (outstr, foutnew)  ;  
    }
    fclose (foutnew) ;
    }
    
    // then write file for those pixels > -.8 but not a new hotfrac alert
    num_alerts_new = alinds.size() ;
    if (num_alerts_new > 0) {
    foutnew = fopen (alertfile, "a") ;
    for (i=0; i< num_alerts_new; i++) {
        
        ind = alinds.at(i) ;
        //zval = (alerts[ind] - basevals_full_new[6*i+4])/basevals_full_new[6*i+5] ;
        zval = (alerts[ind] - glob_nti[ind*2])/glob_nti[ind*2+1] ;
        //zval = alerts[ind] / glob_max21[ind] ;
         //if (zval < 1.1) continue ;
//        sprintf (outstr, "%d\t%s\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%3.1f\r\n", 
//                    ind, modid, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
//                    alerts[ind], zval) ;
        sprintf (outstr, "%d\t%d\t%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%7.3f\t%7.3f\t%7.3f\r\n", 
                    ind, jday, utime, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
                    alerts[ind], newalerts[ind], glob_21[ind*2], glob_32[ind*2]) ;
        //sprintf (outstr, "%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%6.3f\r\n", 
        //            ind, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
        //            alerts[ind], glob_nti[2*ind]) ;
        fputs (outstr, foutnew)  ;  
    }
    fclose (foutnew) ;
    }
    
    fclose (fproc) ;
    delete [] glob_nti;
    delete [] glob_32 ;
    delete [] glob_21 ;
} 

/**
 * calc_alert - The entire modis array, not the resampled grid is run through the alert algorithms. At this
 * time there are two algorithms, the first being if -.8 NTI is exceeded and the second being to compare
 * to the global nti. If the global nti mean plus a threshold number of stdev, then the pixel is alerted.
 * A global alert file for each pixel is generated.
 * @param alertfile
 * @param localertfile
 * @param openFlag
 */
void modis_process::calc_alert(char *alertfile, char *localertfile, bool openFlag) {
    char hdrfile [420], outstr[240], nalertfile[420], *tmp;
    bool ingrid ;
    int i, j, ij, ind, ib, ialert, npix,  num_alerts_new, alert_count[5], ixloc, iyloc, npix_grid;
    float std_alert, zval, val21, val22, val32, val6, alval, alert_limit, new_alertlimit ;
    float  *glob_21, *glob_32, *glob_nti, global_mean, global_std ;
    float z_lat, z_lon, z_y, z_x, xloc, yloc ;
    npix = ns * nl ;
    vector <int> new_alinds ;
    vector <float> std_alvals ;
    vector <float> new_alvals ;
    // temporarily get the local alert file as well
    //strcpy (nalertfile, alertfile) ;
    //tmp = strstr (nalertfile, "_alert.txt") ;
    //strncpy (tmp, "_localert.txt\0",14) ;
    
    npix_grid = nx_grid * ny_grid ;
    
    glob_nti = new float [nl * ns * 2] ;
    glob_21 = new float [nl * ns * 2] ;
    glob_32 = new float [nl * ns * 2] ;
    alinds.clear();
    std_alvals.clear() ;
    new_alvals.clear() ;
    alertFlag = true ;
    for (ialert=0; ialert<4; ialert++) alert_count[ialert] = 0 ;
    alert_limit = night_limit ;
    if (geom->dayflag) {
        alert_limit = day_limit ;
    }
    
    // we get the global nti for each line of data, calc the alerts in that...
    cout << "Getting nti values for the granule" << endl ;
    if (openFlag == true) {
        int hg = 1 ;
        hg * 2 ;
    }
    sp->getAllValues (mon, geom->aq_terra_flag, &geom->latarr[0], &geom->lonarr[0], ns*nl,glob_21, glob_32, glob_nti, openFlag) ;
    cout << "Alerting begun" << endl ;
    
    
    char ntmp[20] ;
    int jday ;
    int utime ;
    
    strncpy (ntmp, &modid[5],3) ;
    ntmp[3]='\0';
    jday = atoi (ntmp) ;
    strncpy (ntmp, &modid[9],4) ;
    ntmp[4] ='\0' ;
    utime = atoi (ntmp) ;
    
    FILE *fproc = fopen (proclog, "a") ;
    for (i=0; i<5; i++)alert_count[i] = 0 ;
    for (i = 0; i < nl; i++) {


        // a pixel exceeding the glob_nti enhanced value is alerted in the new alert algorithm
        for (j = 0; j < ns; j++) {
            
            ij = i * ns + j;
            
            // use mean based upon the global nti values 10yr mean and avg
            global_mean = glob_nti[2*ij] ;
            global_std = glob_nti[2*ij+1] ;
            new_alertlimit = global_mean + 3. * global_std ;
            // actually maxnti
            //new_alertlimit = glob_max21[ij] ;
            if (geom->dayflag){
                if (geom->solzen[ij]>90.) return ;
                
            }
            else {
                if (geom->solzen[ij] < 90) return ;
            }

            //global_mean = glob_21[2*ij] ;
            //global_std = glob_21[2*ij+1] ;
            //new_alertlimit = global_mean + 3. * global_std ;

            // band 6,21,22, 32
            //val6 = therm->refdata_cal[ij];
            val21 = therm->raddata_cal[ij];
            val22 = therm->raddata_cal[npix + ij];
            val32 = therm->raddata_cal[2 * npix + ij];
            
            // check if val22 is saturated, then use val21
            // check if val32 is less than the mean value for 32
            //if (val32 < glob_32[ij*2]) continue ;
            //if either 21 | 22 or 32 < 0. continue 
            if (val22 > 2.0 || val22 < 0.001)
                val22 = val21;
            // also check to make sure 21 and 32 are good 
            if (val21 > 150. || val32 > 40.) {
                *(alerts + ij) = -2.;
                continue;
            }
            if (val22 <= 0. || val32 <= 0. ) {
                *(alerts + ij) = -2.;
                continue ;
            }
            //alval = (val22 - val32) / (val22 + val32);
            std_alert = (val22 - val32) / (val22 + val32);
            alval = std_alert;
            z_lon = geom->lonarr[i * ns + j];
            z_lat = geom->latarr[i * ns + j];
            xloc = (z_lon - startlon) / gridspace;
            yloc = (startlat - z_lat) / gridspace;
            ixloc = int(xloc + 0.5);
            iyloc = int(yloc + 0.5);
            ingrid = true;
            if (ixloc < 0 || iyloc < 0) ingrid = false;
            if (iyloc >= ny_grid || ixloc >= nx_grid) ingrid = false ;
                
            if (alval > new_alertlimit) {
                
                new_alinds.push_back(ij) ;
                
                zval = (alval - global_mean) / global_std ; 
                //zval = (val22 - global_mean) / global_std ;
                if (zval >= 3.0) {
                    alert_count[0]++ ;
                    if (ingrid)
                    histoarr[iyloc*nx_grid+ixloc]+=1 ;
                }
                if (zval >= 4) {
                    alert_count[1]++ ;
                    if (ingrid)
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid]+=1 ;
                }
                    
                if (zval >= 5.) {
                    alert_count[2]++ ;
                    if (ingrid)
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid*2]+=1 ;
                }
                if (zval >= 7.5) {
                    alert_count[3]++ ;
                    if (ingrid)
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid*3]+=1 ;
                }
                if (zval >= 10) {
                    alert_count[4]++ ;
                    if (ingrid)
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid*4]+=1 ;
                }
                
             }
            
            if (std_alert > alert_limit){
                alinds.push_back(ij);
                std_alvals.push_back ( std_alert) ;
                new_alvals.push_back (zval) ;
                if (ingrid) 
                    histoarr[iyloc*nx_grid+ixloc+ npix_grid*5]+=1 ;
            }
            *(alerts + ij) = alval ;
                

        }
    }
    num_alerts_full = alinds.size() ;
    //if (basevals_full) delete[] basevals_full ;
    //basevals_full = new float [6 * num_alerts_full] ;
    cout <<"Number of alerts is modis scene is "  << num_alerts_full << endl ;
    //sp->getCorrespondingValues (alinds, mon, geom->aq_terra_flag, geom->latarr, geom->lonarr, basevals_full, openFlag) ;
    
    if (num_alerts_full > 0) {
    FILE *fout = fopen (alertfile, "a") ;
    for (i=0; i< num_alerts_full; i++) {
        
        ind = alinds.at(i) ;
        //zval = (alerts[ind] - glob_nti[ind*2])/glob_nti[ind*2+1] ;
            sprintf (outstr, "%d\t%d\t%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%6.3f\r\n", 
                    ind, jday, utime, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind], std_alvals[i], new_alvals[i]) ;
                    //std_alvals[i], glob_nti[ind*2], glob_nti[ind*2+1]) ;
        fputs (outstr, fout)  ;  
    }
    fclose (fout) ;
    }
            
    
    num_alerts_new = new_alinds.size() ;
    
//    if (basevals_full_new) delete[] basevals_full_new ;
//    basevals_full_new = new float [6 * num_alerts_new] ;
    
    cout <<"Number of new alerts is modis scene is "  << num_alerts_new << endl ;
    
    //sp->getCorrespondingValues (new_alinds, mon, geom->aq_terra_flag, geom->latarr, geom->lonarr, basevals_full_new) ;
    
    if (num_alerts_new > 0) {
    FILE *foutnew = fopen (localertfile, "a") ;
    for (i=0; i< num_alerts_new; i++) {
        
        ind = new_alinds.at(i) ;
        //zval = (alerts[ind] - basevals_full_new[6*i+4])/basevals_full_new[6*i+5] ;
        zval = (alerts[ind] - glob_nti[ind*2])/glob_nti[ind*2+1] ;
        //zval = alerts[ind] / glob_max21[ind] ;
         //if (zval < 1.1) continue ;
//        sprintf (outstr, "%d\t%s\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%3.1f\r\n", 
//                    ind, modid, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
//                    alerts[ind], zval) ;
        sprintf (outstr, "%d\t%d\t%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%3.1f\t%7.3f\t%7.3f\t%7.3f\r\n", 
                    ind, jday, utime, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
                    alerts[ind], zval, glob_nti[ind*2], glob_21[ind*2], glob_32[ind*2]) ;
        //sprintf (outstr, "%d\t%7.3f\t%8.3f\t%6.2f\t%6.2f\t%6.3f\t%6.3f\r\n", 
        //            ind, geom->latarr[ind], geom->lonarr[ind], therm->raddata_cal[ind], therm->raddata_cal[2*npix+ind],
        //            alerts[ind], glob_nti[2*ind]) ;
        fputs (outstr, foutnew)  ;  
    }
    fclose (foutnew) ;
    }
    
   
    sprintf (outstr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\r\n", jday, utime, num_alerts_full, alert_count[0], alert_count[1], alert_count[2], 
            alert_count[3], alert_count[4]) ;
    fputs (outstr, fproc) ;
    
    fclose (fproc) ;
    alinds.clear() ;
    delete [] glob_nti;
    delete [] glob_32 ;
    delete [] glob_21 ;
    
}


   
      
// write out each alert to an appended file with the following columns
// mod21 filename
// year
// month
// day
// sol zenith
// b21
// b22
// b32
// alert value

void modis_process::write_alert_textfile(char *outfile) {
    FILE *fout = fopen(outfile, "a+");
    int nalerts = alinds.size();

}

void modis_process::write_output(char *outfile) {
    char hdrfile [420];
    int i, npix;
    npix = nx_grid * ny_grid;
    cout << "Number of samples is : " << nx_grid << endl;
    cout << "Number of lines is : " << ny_grid << endl;

    strcpy(hdrfile, outfile);
    strcat(hdrfile, ".hdr");
    write_header(hdrfile, 6);
    cout << "Output file is : " << outfile << endl ;

    FILE *fout = fopen(outfile, "w");
    if (fout == NULL) {
        cout << "Could not open " << outfile << endl;
        return;
    }
    
    
    fwrite((char *) bandsarr, 4, 7 * npix, fout);
    
    
    // band21,22,32,alert, mean stdev
    //fwrite((char *) this->glob_mnstdev, 4, 6 * npix, fout) ;
    
    
   // fwrite ((char *) this->bandsarr, 4, 8* npix, fout) ;
    
    fclose(fout);

}

void modis_process::write_header(char *outfile, int nbands) {

    FILE *hdrout = fopen(outfile, "w");
    char bnames [1200];
    strcpy(bnames, "band names = {\nModisBand21,ModisBand22,ModisBand32,Alert,Alert_z, Alert_rat}\n") ;
            //"10yr_2122,10yr_2221_stdv,10yr_32,10yr_32_stdv, 10yr_nti, 10yr_nti_stdv}\n");

    fprintf(hdrout, "ENVI\ndescription = {\nMOD021KM - MOD03  }\n");
    fprintf(hdrout, "samples    = %5d\n", nx_grid);
    fprintf(hdrout, "lines      = %5d\n", ny_grid);
    fprintf(hdrout, "bands      = %3d\n", nbands);
    fprintf(hdrout, bnames);
    fprintf(hdrout, "header offset = 0 \n");
    fprintf(hdrout, "file type = ENVI Standard \n");
    fprintf(hdrout, "data type = 4 \n");
    fprintf(hdrout, "interleave = bsq \n");
    fprintf(hdrout, "sensor type = MODIS \n");
    fprintf(hdrout, "byte order = 0\n");
    fprintf(hdrout, "map info = {Geographic Lat/Lon, 1.0000, 1.0000, %12.8f, %12.8f,",
            startlon, startlat);
    fprintf(hdrout, "%6.4f, %6.4f, WGS-84, units=Degrees}\n", gridspace, gridspace);

    fclose(hdrout);
}

float modis_process::bb_temptorad (float wave, float temp) {
    double rad, den ;
    //float tempK = temp + 273.15 ;
    float tempK = temp ;
    rad = 1.191066e8 / pow(wave,5.) ;
    den = (exp(1.4388E4/(wave *tempK))-1.) ;
    return float(rad/den) ;
}

float modis_process::bb_radtotemp(float wave, float rad) {

    float wave_m, l_m;
    double h, k, c, c1, c2, val0, Temp;

    // convert to m
    l_m = rad * 1.E6;
    wave_m = wave * 1.E-6;
    h = 6.62606755E-34;
    k = 1.380658E-23;
    c = 2.9979246E8;
    c1 = 2. * h * c * c;
    c2 = h * c / k;

    val0 = c2 / (wave_m * log(c1 / (l_m * pow(wave_m, 5)) + 1.));

    return float(val0);
}

void modis_process::zscore () {
    int i, npix ;
    npix = nx_grid * ny_grid ;
    float alval ;
    for (i=0; i<npix; i++) {
        
        alval = (bandsarr[3*npix+i]-glob_mnstdev[4*npix+i]) / glob_mnstdev[5*npix+i] ;
        // zscore
        bandsarr[4*npix+i] = alval ;
        // alert ratio
        bandsarr[5*npix+i] = bandsarr[3*npix+i] / glob_mnstdev[4*npix+i] ;
    }
    
    
}




void modis_process::extract_from_baseline_file () {
    int i ;
    
    
    //sp->set_projection_parameters ()
    sp->fillGrid (this->mon, this->geom->aq_terra_flag, startlat, startlon, gridspace, nx_grid, ny_grid, glob_mnstdev) ;
    //for (i=0; i<nx_grid*ny_grid; i++) {
    
}



/*
https://cimss.ssec.wisc.edu/dbs/China2011/Day2/Lectures/MODIS_DB_Land_Surface_Temperature_reference.pdf
float modis_process::bb_radtotemp (float wave, float rad) {

        float wave_m, l_m ;
        double h, k, c, K1, K2, val0, Temp ;

        // wave in microns 
        l_m = rad   ;
        wave_m = wave  ;
        K2 = 1.43883E4 / wave;
        K1 = 1.19107E8 ;
        val0 = K1 / (l_m*pow(wave,5)) + 1. ;
        Temp = (K2 / log (val0))  ;
        return float(Temp) ;
}


 */
