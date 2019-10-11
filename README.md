# modis_alert_dev
Suite of programs  to arrive at a new MODIS alert algorithm. Generally, the algorithm will calculate the global limit, whether it 
be a limit for band 21/22 (4 micron) or for a NTI value (b221-b32)/(b221+b32). 
Each incoming MODIS pixel will search the global map to determine if it exceeds the limit, and hence is identified as an alerted 
pixel.

# contents and operation
Numerous .h and .cpp files, with a Makefile which facilitates the building of the various programs. Only requirement is the hdf-4.2.13 library. 
