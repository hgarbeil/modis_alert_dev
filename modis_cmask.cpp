#include "modis_cmask.h"

modis_cmask::modis_cmask (const char *ifile) {

	ns = 1354 ;
	nl = 2030 ;
	cmask = new unsigned char [ns *nl] ;
	int status = load_cmask(ifile) ;

}
modis_cmask::~modis_cmask () {
	delete [] cmask ;
}



int modis_cmask::load_cmask(const char *ifile) {
	int i, sd_id,n_datasets, n_fileattrs; 
	int32 sds_id, rank, dim_sizes[5], attributes, num_type, attnum ;
	char name [240], attr_name[40] ;
	int32 stride[3]={1,1,1} ;
	int32 edge[3]={1,2030,1354} ;
	int32 start[3]={0,0,0} ;

	FILE *fin = fopen (ifile, "r") ;
	if (fin==NULL) {
		cout <<"Could not open "<< ifile<< endl ;
		return (-1) ;
	}
	fclose (fin) ;
	sd_id = SDstart (ifile, DFACC_RDONLY) ;
	if (sd_id <= 0){
		cout << "Problem with : " << ifile << endl ;
		return -1 ;
	}
	cout << "opened " << ifile << " with id : " << sd_id << endl ;


	SDfileinfo (sd_id, &n_datasets, &n_fileattrs) ;
	//cout << "num datasets " << n_datasets << " file atts : " << 
	//	n_fileattrs << endl ;
	sds_id = SDselect (sd_id, 8) ;
	SDgetinfo (sds_id, name, &rank, dim_sizes, &num_type, &attributes) ;
	SDreaddata (sds_id, start, stride, edge, cmask) ;
	SDendaccess(sds_id) ;

	SDend (sd_id) ;
	FILE *fout = fopen ("cloudmask", "w") ;
	fwrite (cmask, 1, 1354 * 2030, fout) ;
	fclose (fout) ;


	return (sd_id) ;
}

void modis_cmask::make_map (unsigned char *map ){

	int i ;
	unsigned char uval, outval, tval;

	for (i=0; i<ns *nl;i++) {
		outval = 0 ;
		if (i==(1354*521+750)) {
			int hg = 0 ;
			hg *=1 ;
		}
		uval = cmask [i] ;

		tval = (uval>>1)&1 ;
		outval = outval + tval ;
		tval = (uval>>2)&1 ;
		outval = outval + tval * 10 ;
		map[i]=outval ;
	}

}	
/*
		
int main (int argc, char *argv[]) {

	modis_cmask *mc = new modis_cmask ("/local/worldbase/MOD35_L2/2001/001/MOD35_L2.A2001001.2100.061.2017219155751.hdf") ;

	unsigned char *map=new unsigned char [1354 * 2030] ;
	mc->make_map(map) ;


	delete mc ;
	delete [] map ;
}
*/
