#include "alert.h"
#include "math.h"
#include <stdio.h>


alert::alert() {

	// standard size of modis granule, we will have a method to switch to viirs sizes 
	ns = 1354 ;
	nl = 2030 ;
	viirsflag = false  ;
	npix = ns * nl ;

	nti_thresh = -.80 ;
	badpix = 0l ;
	b21max = 100. ;
	b22max = 2.2 ;

}
void alert::set_max (float b21v, float b22v) {
	b21max = b21v ;
	b22max = b22v ;
}

void alert::set_badpix (unsigned char *b) {
	badpix = b ;
}

void alert::set_viirs(bool vflag) {
	if (vflag) {
		nl = 3000 ;
		ns = 3224 ;
		viirsflag = true ;
	}
	else {
		ns = 1354 ;
		nl = 2030 ;
		viirsflag = false  ;
	}
	npix = ns * nl ;
}


// based on band 32 avg or current, get band 22 radiance but add 10degrees to TB4micron -> L4
int alert::calc_varnti (float *b21, float *b22, float *b32, float *exval,  float *mnstdv, float *mnstdv_32, vector <int> &al_bbind, vector <float> &al_bbvec) {
        
        int i, count=0 ;
	
        float   b22high, b22pred, meanval, stdval,b21val, b32val, b22val, b221val, exceed_val, hot_rad_frac, temp12 ;
	float b22_alert, b22_exceed, temp4, ntival, ntilimit ;
	float minfrac, fraclevels[] = {.00015, .000175, .0002, .00025, .0003} ;



	for (i=0; i<npix; i++) {
		exval[i] =0. ;
		if (badpix[i]) {
                    //cout << "bad pix hit" <<endl ;
                    continue ;
		}
		b32val = b32[i] ;
        b22val = b22[i] ;

		// get a high edge for band 22 for this location
		meanval = mnstdv[i*2] ;
		stdval = mnstdv[i*2+1] ;
		b22high = meanval + 2. * stdval ;
        if (b22val < b22high) continue ;
		//b22mean = meanval ;
		// get the greater of the 10 yr avg or the img value
		// ** temporary turn off

        // get the warmer of the two band32 temperatures, this or global
		meanval = mnstdv_32[i*2]+0.*mnstdv_32[i*2+1] ;
		if (b32val < meanval) b32val = meanval ;
		//
		// from the band32 value get the temperature
		//
		// from the band32 value get the temperature
		temp12 = bb_radtotemp (12.0, b32val) ;
		ntilimit = 3E-5 * temp12 * temp12 -.0114*temp12 + 0.22 - 0.3 ;

		// then the radiance from temp12 + a temp difference of X degrees
		//b22pred = bb_temptorad (3.95, temp12+5.) ; 
		
		b21val = b21[i] ;
		b221val = b22val ;
		// bad vals
		if (b22val <= 0  && b21val <=0) continue ;
		// saturated for both 4 micron bands
		if (b22val > b22max && b21val > b21max) continue ;
		
		if (b22val >  b22max && b21val >0){  
			b221val = b21val ;
		}
        // if (b221val < b22pred) continue ;
        //temp4 = bb_radtotemp (3.95, b221val) ;
		// calc the nti for this pixel
		ntival = (b221val - b32val)/(b21val + b32val) ;
		if (ntival <= ntilimit) continue ;
		exval[i] = ntival - ntilimit ;
		al_bbind.push_back(i) ;
		al_bbvec.push_back(ntival) ;
		// save exval as the background value
	}

	return al_bbind.size() ;
}    

		


// based on band 32 avg or current, get band 22 radiance but add 10degrees to TB4micron -> L4
int alert::calc_bb1 (float *b21, float *b22, float *b32, float *exval,  float *mnstdv, float *mnstdv_32, vector <int> &al_bbind, vector <float> &al_bbvec) {
        
        int i, count=0 ;
	
        float   b22high, b22pred, meanval, stdval,b21val, b32val, b22val, b221val, exceed_val, hot_rad_frac, temp12 ;
	float b22_alert, b22_exceed, temp4  ;
	float minfrac, fraclevels[] = {.00015, .000175, .0002, .00025, .0003} ;



	for (i=0; i<npix; i++) {
		exval[i] =0. ;
		if (badpix[i]) {
                    //cout << "bad pix hit" <<endl ;
                    continue ;
		}
		b32val = b32[i] ;
        b22val = b22[i] ;

		// get a high edge for band 22 for this location
		meanval = mnstdv[i*2] ;
		stdval = mnstdv[i*2+1] ;
		b22high = meanval + 2.5 * stdval ;
        if (b22val < b22high) continue ;
		//b22mean = meanval ;
		// get the greater of the 10 yr avg or the img value
		// ** temporary turn off

        // get the warmer of the two band32 temperatures
		meanval = mnstdv_32[i*2]+2.*mnstdv_32[i*2+1] ;
		if (b32val < meanval) b32val = meanval ;
		//
		// from the band32 value get the temperature
		temp12 = bb_radtotemp (12.0, b32val) ;
		//ntilimit = 3E-5 * temp12 * temp12 -.0114*temp12 + 0.22 ;

		// then the radiance from temp12 + a temp difference of X degrees
		//b22pred = bb_temptorad (3.95, temp12+5.) ; 
		
		b21val = b21[i] ;
		b221val = b22val ;
		//alval[i] = -5.;
		
		// bad vals
		if (b22val <= 0  && b21val <=0) continue ;
		// saturated for both 4 micron bands
		if (b22val > b22max && b21val > b21max) continue ;
		
		if (b22val >  b22max && b21val >0){  
			b221val = b21val ;
		}
        // if (b221val < b22pred) continue ;
        //temp4 = bb_radtotemp (3.95, b221val) ;

		exval[i] = temp4 - temp12 ;
		if (exval[i] < 10.) continue ;
		
		//if (b221val < b22pred) continue ;
		/*
		exceed_val = b221val - b22pred ;
		hot_rad_frac = exceed_val / b22add ;
		*/
		//if (hot_rad_frac > minfrac) {
		
		al_bbind.push_back(i) ;
		al_bbvec.push_back(exval[i]) ;
		// save exval as the background value
		exval[i] = temp12 ;

	}

	return al_bbind.size() ;
}    



// use the band 32 value to predict the band 21 value based upon planck function, then see if the actual 21 exceeds that 
int alert::calc_bb (float *b21, float *b22, float *b32, float *exval,  float *mnstdv, float *mnstdv_32, vector <int> &al_bbind, vector <float> &al_bbvec) {
	int i, count=0 ;
	float b22mean, b22high, stdval, b22add, b22pred, meanval, b21val, b32val, b22val, b221val, exceed_val, hot_rad_frac, tempbb ;

	float b22_alert, b22_exceed ;
	float minfrac, fraclevels[] = {.00015, .000175, .0002, .00025, .0003} ;
	b22add  = this->bb_temptorad(3.95,500+273.15) ;


	minfrac = .00005 ;



	for (i=0; i<npix; i++) {
		exval[i] =0. ;
		if (badpix[i]) {
		//cout << "bad pix hit" <<endl ;
		continue ;
		}
		b32val = b32[i] ;

		// get a high edge for band 22 for this location
		meanval = mnstdv[i*2] ;
		stdval = mnstdv[i*2+1] ;
		b22high = meanval + 3.0 * stdval ;
		b22mean = meanval ;
		// get the greater of the 10 yr avg or the img value
		// ** temporary turn off

		meanval = mnstdv_32[i*2] ;
		if (b32val < meanval) b32val = meanval ;
		//
		// from the band32 value get the temperature
		tempbb = bb_radtotemp (12.0, b32val) ;
		// then the radiance
		b22pred = bb_temptorad (3.95, tempbb+0.) ; 
		if (b22pred < b22mean) b22pred = b22mean ;
		
		b22val = b22[i] ;
		b21val = b21[i] ;
		b221val = b22val ;
		//alval[i] = -5.;
		if (b32val < 0) continue ;
		// bad vals
		if (b22val <= 0  && b21val <=0) continue ;
		// saturated
		if (b22val > b22max && b21val > b21max) continue ;
		
		if (b22val >  b22max && b21val >0){  
			b221val = b21val ;
		}
		b22_exceed = b221val - b22pred ;
		b22_exceed /= b22add ;
		
		//if (b221val < b22pred) continue ;
		/*
		exceed_val = b221val - b22pred ;
		hot_rad_frac = exceed_val / b22add ;
		*/
		//if (hot_rad_frac > minfrac) {
		if (b221val > b22high && b22_exceed > .0001 ) 
		{
			al_bbind.push_back(i) ;
			al_bbvec.push_back(b22_exceed) ;
			exval[i] = b22_exceed ;
		}	

	}

	return al_bbind.size() ;
}

int alert::calc_nti (float *b21, float *b22, float *b32, float *alval, vector <int> &alind, vector <float> &alvec) {
	int i, count=0 ;
	float b21val, b32val, b22val, b221val ;


	for (i=0; i<npix; i++) {
		if (badpix[i]) {
		//cout << "bad pix hit" <<endl ;
		continue ;
		}
		b32val = b32[i] ;
		b22val = b22[i] ;
		b21val = b21[i] ;
		b221val = b22val ;
		alval[i] = -5.;
		if (b32val < 0) continue ;
		if (b22val <= 0  && b21val <=0) continue ;
		if (b22val > b22max && b21val > b21max) continue ;
		
		if (b22val >  b22max && b21val >0){  
			b221val = b21val ;
		}
		alval[i] = (b221val - b32val) / (b221val + b32val) ;
		if (alval[i] > nti_thresh) {
			alind.push_back(i) ;
			alvec.push_back(alval[i]) ;
		}
	}

	return alind.size() ;
}


/***
 * calc_nstdv - uses the global mean and std dev for band 221, checks the b22 and b21 array to arrive at the 4 micron radiance
 * and then checks to see if this value exceeds the mean plus 3 standard deviations. 
 * If so, this is an alerted pixex and adds the number of std devs to the alind and alvec vectors
 * Note that the al_std is the zval for each pixel, not just vecteor of alerted pixels.
 ***/
int alert::calc_nstdv (float *b21, float *b22, float *mnstdv, float *al_std, vector <int> &alind, vector <float> &alvec) {

	int i, count=0 ;
	float zval, b21val, b32val, b22val, b221val, meanval, stdval ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		b21val = b21[i] ;
		b22val = b22[i] ;
		al_std[i] = 0. ;
		meanval = mnstdv[i*2] ;
		stdval = mnstdv[i*2+1] ;
		if (stdval <= 0.000001) continue ;
		b221val = b22val ;
                                    // check for bad data
		if (b21val < 0 && b22val <0) continue ;
                                    // check for saturation
		if (b22val > b22max && b21val > b21max) continue ;
		if (b22val > 2.2) {
			b221val = b21val ;
		}
		zval =( b221val - meanval) / stdval ;
		al_std[i] = zval ;
		if (zval >3.0) {
			alind.push_back(i) ;
			alvec.push_back (zval) ;
			//cout << alind.size() << " " << zval << endl ;
		}

	}
	return alind.size() ;
}
/***
* calc_max 
* method to go through the modis b22 values to check if they exceed the global values found in max2232. If so, the 
* values for that pixel are appended to the relevant vectors, alind and alvec. 
***/

int alert::calc_max (float *b21, float *b22, float *b32, float *max2232, float *al_std, vector <int> &alind, vector <float> &alvec) {

	int i, count=0 ;
	float zval, b21val, b32val, b22val, b221val, maxval ;
	float T32, T22 ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		b21val = b21[i] ;
		b22val = b22[i] ;
		b32val = b32[i] ;
		al_std[i] = 0. ;
		maxval = max2232[i] ;
		b221val = b22val ;
                                    // check for bad data
		if (b21val < 0 && b22val <0) continue ;
        // check for saturation
		if (b22val > b22max && b21val > b21max) continue ;
		// only use 21 if 22 is saturated
		if (b22val > 2.0) {
			b221val = b21val ;
		}
		T22 = bb_radtotemp (4, b22val) ;
		T32 = bb_radtotemp (12, b32val) ;
		if (T22-T32 <15. && zval < 1.15) continue ;
		zval =( b221val / maxval) ;
		al_std[i] = zval ;
		if (zval >.95) {
			alind.push_back(i) ;
			alvec.push_back (zval) ;
			//cout << alind.size() << " " << zval << endl ;
		}

	}
	return alind.size() ;

}


int alert::calc_nti_stdv (float *b21, float *b22, float *b32, float *mnstdv, float *al_std, vector <int> &alind, vector <float> &alvec) {

	int i, count=0 ;
	float ntival, zval, b21val, b32val, b22val, b221val, meanval, stdval ;

	for (i=0; i<npix; i++) {
		if (badpix[i]) continue ;
		b21val = b21[i] ;
		b22val = b22[i] ;
                b32val = b32[i] ;
		al_std[i] = 0. ;
		meanval = mnstdv[i*2] ;
		stdval = mnstdv[i*2+1] ;
		b221val = b22val ;
		if (b21val < 0 && b22val <0) continue ;
		if (b22val > b22max && b21val > b21max) continue ;
		if (b22val > 2.0) {
			b221val = b21val ;
		}
                ntival = (b221val - b32val) / (b221val + b32val) ;
		zval =( ntival - meanval) / stdval ;
		al_std[i] = zval ;
		if (zval > 3.0) {
			alind.push_back(i) ;
			alvec.push_back (zval) ;
			//cout << alind.size() << " " << zval << endl ;
		}

	}
	return alind.size() ;
}

// the alert grid file has 5 bands, 1st band is mean radiance, 2 is
// old alerthis, 3 2stdev, 4 2.6 stdv, 5 3sdev
void alert::write_envi_header (const char *outfile, int ns, int nl, float start_lat, float start_lon, float gspace) {
	FILE *fhdr = fopen (outfile, "w") ;
	fprintf(fhdr, "Envi = \r\ndescription = {\r\nMODIS Alert Grid Results.}\r\n") ;
	fprintf(fhdr, "samples = %d\r\n", ns) ;
	fprintf(fhdr, "lines = %d\r\n", nl) ;
	fprintf(fhdr, "bands = 6\r\n") ;
	fprintf(fhdr, "file type = ENVI Standard\r\n") ;
	fprintf(fhdr, "data type = 12\r\n") ;
	fprintf(fhdr, "interleave = bsq\r\n") ;
	fprintf(fhdr, "band names ={Radiance,RadianceStdv,NTI Hits,3_SDev, 4_SDev, 5_SDev}\r\n") ;
	fprintf (fhdr, "map info ={Geographic Lat/Lon, 1., 1., %f, %f, %f, %f, WGS-84, units=Degrees",
		start_lon, start_lat, gspace, gspace) ;

	fclose (fhdr) ;
}

// the alert grid file has 5 bands, 1st band is mean radiance, 2 is
// old alerthis, 3 2stdev, 4 2.6 stdv, 5 3sdev
void alert::write_envi_header (const char *outfile, int ns, int nl, float start_lat, float start_lon, float gspace, int dtype) {
	FILE *fhdr = fopen (outfile, "w") ;
	fprintf(fhdr, "Envi = \r\ndescription = {\r\nMODIS Alert Grid Results.}\r\n") ;
	fprintf(fhdr, "samples = %d\r\n", ns) ;
	fprintf(fhdr, "lines = %d\r\n", nl) ;
	fprintf(fhdr, "bands = 3\r\n") ;
	fprintf(fhdr, "file type = ENVI Standard\r\n") ;
	fprintf(fhdr, "data type = %d\r\n", dtype) ;
	fprintf(fhdr, "interleave = bsq\r\n") ;
	fprintf(fhdr, "band names ={3_SDev, 4_SDev, 5_SDev}\r\n") ;
	fprintf (fhdr, "map info ={Geographic Lat/Lon, 1., 1., %f, %f, %f, %f, WGS-84, units=Degrees",
		start_lon, start_lat, gspace, gspace) ;


	fclose (fhdr) ;

}

float alert::bb_radtotemp(float wave, float rad) {

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

float alert::bb_temptorad (float wave, float temp) {
     double rad, den ;

     float tempK = temp ;
     rad = 1.191066e8 / pow(wave,5.) ;
     den = (exp(1.4388E4/(wave *tempK))-1.) ;
     return float(rad/den) ;
}
													
		
	
