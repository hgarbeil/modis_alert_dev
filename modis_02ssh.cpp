
#include "modis_02ssh.h"
#define DtoR M_PI / 180. 

modis_02ssh::modis_02ssh () {
    e_radius = 6371007.181 ;
    start_x = -1. * M_PI * e_radius ;
    start_y = M_PI * e_radius * 0.5 ;
    space = 4633.1271655 ;
    
    
    
}

void modis_02ssh::get_mnstdv (float *lon, float *lat, float *mnstdv, int npts) {
    int i ;
    float stdval, mnval ;
    
    for (i=0; i<npts; i++){
        latval = lat[i] * DtoR ;
        lonval = lon[i] *DtoR ;
        
        
        
        
        
    }
    
    
}
