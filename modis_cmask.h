#ifndef mcmask
#define mcmask

#include <fstream>
#include <iostream>
#include <stdio.h>
#include "mfhdf.h"


using namespace std ;


class modis_cmask {
	public :
	modis_cmask (const char *ifile);
	~modis_cmask () ;
	int load_cmask (const char *) ;
	void make_map (unsigned char *) ;
	char inname [420] ;
	unsigned char *cmask;
	int ns, nl ;

};

#endif
