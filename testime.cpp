#include <time.h>
#include <stdio.h>


long converttime_unix (int year, int month, int day, int hr, int min, int sec) {
	struct tm t ;
	time_t t_unix ;

	t.tm_year = year - 1900 ;
	t.tm_mon = month - 1 ;
	t.tm_mday = day ;
	t.tm_hour = hr ;
	t.tm_min = min ;
	t.tm_sec = int(sec) ;

	t_unix = mktime (&t) ;
	//printf ("Unixtime is %ld\r\n", (long)t_unix-36000) ;
	return (long)(t_unix-36000) ;
}
