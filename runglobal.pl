#!/usr/bin/perl -s

for (my $day=200;$day<250;$day++) {

$dir=sprintf("/hotspot3/data2/maxalerts/2014/%03d",$day) ;
mkdir $dir;
$dir1=sprintf("/hotspot3/data2/alerts_exclusive/2014/%03d",$day) ;
mkdir $dir1;

$cmd = sprintf("modis_gmax globfiles/day_%03d.txt results_global/terra_night_%03d 3 %03d", $day, $day, $day) ;
system($cmd) ;
}
