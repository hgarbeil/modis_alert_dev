//#include "modis_hdf"
#include <sys/stat.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <glob.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "modis_cmask.h"
using std::vector;
using namespace std;

vector<string> globVector(const string& pattern) ;
void get_line_samp (const char *, int *, int *) ;

int main (int argc, char *argv[]) {
	bool m35_notfound ;
	int i, iday, lines, index, yloc, xloc ;
	char indir[240], inrec[420], in35[420], outdir0[420] ;
	unsigned char *cmask ;
	string fname, cfile ="", outfilenm ;
	string line, matchstring, outdir ;
	vector<string> mfiles ;
	vector<string> m35files ;
	ifstream alfile ;
	ofstream outfile ;
	modis_cmask *mc=NULL ;
	stringstream numstr;

	//strcpy (indir, *++argv) ;
	// get the list of alert files
	//
	cmask = new unsigned char [1354 * 2030] ;
	for (iday=0; iday<366;iday++) {
	sprintf (indir, "/home/modis/done/done2002/%03d/M*dat",iday) ;
	sprintf (inrec, "/local/worldbase/MOD35_L2/2002/%03d/M*hdf",iday) ;
	sprintf (outdir0, "/user3/hg/modvolc/2002/%03d/", iday) ;
	mkdir (outdir0, 0755) ;

	
	//mfiles = globVector ("/home/modis/done/done2001/001/M*dat") ;
	mfiles.clear() ;
	mfiles = globVector (indir) ;
	//m35files = globVector ("/local/worldbase/MOD35_L2/2001/001/M*hdf") ;
	m35files.clear() ;
	m35files = globVector (inrec) ;
	outdir =  outdir0 ;
	cout << "number of files is : " << mfiles.size() << endl ;
	//mfiles.clear() ;
	//mfiles.push_back ("/home/modis/done/done2001/001/MODVOLC.A2001001.0335.000.000000.dat") ;

	for (i=0; i<mfiles.size() ; i++) {
		index = mfiles[i].find ("A20") ;
		fname = mfiles[i].substr (index-8, string::npos) ;
		outfilenm = outdir+fname ;
		cout << "file name is : " << outfilenm << endl ;
		
		alfile.open (mfiles[i].c_str(), ifstream::in) ;
		outfile.open (outfilenm.c_str(), ofstream::out) ;
		for (lines=0; getline(alfile,line); lines++) {
			if (lines==0) 
				outfile << line << endl ;
		}
		alfile.close() ;
		cout << "Number of lines is " << lines << endl ;
		if (lines==1) {
			outfile.close() ;
			continue ;
		}
		
		m35_notfound = true ;
		index = mfiles[i].find("A20") ;
		matchstring = mfiles[i].substr (index, 13) ;
		
		for (vector<string>::iterator tfile = m35files.begin() ; 
			tfile!=m35files.end() ; tfile++) {
			cfile = *tfile ;	
			index = cfile.find(matchstring)  ;
			if (index>0) {
				cout << "matching  : " << cfile << endl ;
				m35_notfound = false ;
				break ;
				
			}
		}
		if (!m35_notfound){
			if (mc) delete mc ;
			mc= new modis_cmask (cfile.c_str()) ;
			mc->make_map(cmask) ;
		}
		string test;
		alfile.open (mfiles[i].c_str(), ifstream::in) ;
		getline (alfile,line) ;
		for (lines=0; getline(alfile,line); lines++) {
			get_line_samp (line.c_str(), &yloc, &xloc) ;
			//numstr << setw(2) << setfill('0') << cmask[yloc*1354+xloc] ;
			if (!m35_notfound)
			{
				if (cmask[yloc*1354+xloc]<2) 
					line = line + " 0"+to_string(cmask[yloc*1354+xloc]);
				else 
					line = line + " "+to_string(cmask[yloc*1354+xloc]);
			}
			else
				line = line + " -99";
			outfile << line  << endl ;
			
		}
		alfile.close() ;
		outfile.close() ;
		
	} // end of file

	} // end of days

	delete [] cmask ;
	if (mc) delete mc ;
		
	



}

vector<string> globVector(const string& pattern){
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    	files.push_back(string(glob_result.gl_pathv[i]));
	}
	globfree(&glob_result);
	return files;

}


void get_line_samp (const char *instr, int *line, int *samp) {

	int i ;
	char tempstr[240], *tmp ;
	strcpy (tempstr, instr) ;

	tmp = strtok (tempstr, " ") ;
	for (i=0; i<10; i++) {
		tmp = strtok (NULL, " ") ;

	}
	*line = atoi (tmp) ;
	tmp = strtok (NULL, " ") ;
	*samp = atoi (tmp) ;
	
}


