#ifndef aldef
#define aldef
#include <iostream>
#include <vector>

using namespace std ;

class alert {
public :
	alert () ;
	unsigned char *badpix ;
	void set_max (float, float) ;
	void set_badpix (unsigned char *b) ;
	int calc_nti (float *b21, float *b22, float *b31, float *alval, vector <int>&alind, 
		vector<float>&alvec) ;
	int calc_nstdv (float *b21, float *b22, float *mnstdv, float *alval, vector <int> &alind, vector <float> &alvec) ;
        int calc_nti_stdv (float *b21, float *b22, float *b32, float *mnstdv, float *al_std, vector <int> &alind, vector <float> &alvec)  ;
	int calc_bb(float *b21, float *b22, float *b32, float *exval, float*mnstdv, float *mnstdv_32, std::vector<int>&, std::vector<float>&) ;
	int calc_bb1 (float *b21, float *b22, float *b32, float *exval, float*mnstdv, float *mnstdv_32, std::vector<int>&, std::vector<float>&) ;
	int calc_varnti (float *b21, float *b22, float *b32, float *exval, float*mnstdv, float *mnstdv_32, std::vector<int>&, std::vector<float>&) ;
	void write_envi_header (const char *outfile, int ns, int nl, float start_lat, float start_lon, float gspace) ;
	void write_envi_header (const char *outfile, int ns, int nl, float start_lat, float start_lon, float gspace, int dt) ;
	float bb_radtotemp (float , float) ; 
	float bb_temptorad (float , float) ; 
	void set_viirs (bool) ;	
	int ns, nl, npix ;
	float b21max, b22max ;
	bool viirsflag ;
	float nti_thresh ;

} ;

#endif
