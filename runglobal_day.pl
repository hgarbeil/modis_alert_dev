#!/usr/bin/perl -s

for (my $day=120;$day<180;$day++) {

$dir=sprintf("/hotspot3/data2/maxalerts_day/2014/%03d",$day) ;
mkdir $dir;
$dir1=sprintf("/hotspot3/data2/alerts_exclusive_day/2014/%03d",$day) ;
mkdir $dir1;

$cmd = sprintf("modis_gmax_day globfiles/day_%03d.txt results_global_day/terra_day_%03d 1 %03d", $day, $day, $day) ;
system($cmd) ;
}
