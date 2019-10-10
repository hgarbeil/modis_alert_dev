#include "alert.h"
#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {

	int i ;
	float nti, tempbb, addtemp, radbb12, radbb4, radbb4hot, temp12, temp4 ;

	alert *al = new alert() ;
	for (i=0; i<10; i++) {
		tempbb = 233 + i * 10. ;
		radbb12 = al->bb_temptorad (12., tempbb) ;
		//radbb12 *= .75 ;
		addtemp =  400. ;
		radbb4 = al->bb_temptorad (4., tempbb+10.) ;
		radbb4hot = al->bb_temptorad (4., addtemp) ;
		radbb4 = .9999 * radbb4 + .0001*radbb4hot ;
		nti = (radbb4-radbb12)/ (radbb4+radbb12) ;
		cout << tempbb << "\t"<<   nti << endl ;
		temp12 = al->bb_radtotemp(12,3.93) ;
		temp4 = al->bb_radtotemp(3.95,.277) ;
		nti = (.277 - 3.93)/(.277+3.93) ;
		//cout << temp12 << "\t"<< temp4 << "  " << nti << endl ;

	}
	delete al ;
}

