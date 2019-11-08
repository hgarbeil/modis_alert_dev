CC = g++ -std=c++1y -g

# for linux
hdflibs = -L/hbeta/harold/lhome/external/hdf-4.2.13/lib -lm -lmfhdf -ldf -ljpeg -lz -lsz -lboost_system -lboost_filesystem
hdfinc = -I/hbeta/harold/lhome/external/hdf-4.2.13/include 
# for macos 
#hdflibs = -L/Users/hg/workdir/hdf-4.2.5-mac-intel-x86_64/lib -lm -lmfhdf -ldf -ljpeg -lz -lsz -lboost_system -lboost_filesystem
#hdfinc = -I/Users/hg/workdir/hdf-4.2.5-mac-intel-x86_64/include

OBJS = main.o  modis_hdf.o  surftemp.o modis_process.o
OBJS_alert = main_alert.o  modis_hdf.o  surftemp.o modis_process.o sinuProjection.o
OBJS_alert_max = main_max.o  testime.o modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o
OBJS_histo = main_histo.o  testime.o modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o
OBJS_stats = main_stats.o  modis_hdf.o  testime.o surftemp.o modis_process.o alert.o sinu_1km.o sinuProjection.o
OBJS_stats_sinu = main_stats_sinu.o  modis_hdf.o  testime.o surftemp.o modis_process.o alert.o sinu_1km.o sinuProjection.o
OBJS_maxtile = main_maxtile.o  modis_hdf.o  testime.o surftemp.o modis_process.o alert.o sinuProjection.o sinu_1km.o
OBJS_bb = main_bb.o  modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o testime.o
OBJS_varnti = main_varnti.o  modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o testime.o
OBJS_base = main_base.o  sinu_1km.o surftemp.o
OBJS_tseries = main_tseries.o  modis_hdf.o
OBJS_alert_nti_dev = main_nti.o  testime.o modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o
OBJS_alert_sdev = main_sdev.o  testime.o modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o
OBJS_histo = main_histo.o  testime.o modis_hdf.o  surftemp.o modis_process.o sinu_1km.o sinuProjection.o alert.o

modis_maxtile : ${OBJS_maxtile}
	${CC} ${OBJS_maxtile} ${hdflibs} -o modis_maxtile

modis_stats : ${OBJS_stats}
	${CC} ${OBJS_stats} ${hdflibs} -o modis_stats

modis_stats_sinu : ${OBJS_stats_sinu}
	${CC} ${OBJS_stats_sinu} ${hdflibs} -o modis_stats_sinu

modis_nti_dev : ${OBJS_alert_nti_dev}
	${CC} ${OBJS_alert_nti_dev} ${hdflibs} -o modis_nti_dev

modis_histo : ${OBJS_histo}
	${CC} ${OBJS_histo} ${hdflibs} -o modis_histo

modis_sdev : ${OBJS_alert_sdev}
	${CC} ${OBJS_alert_sdev} ${hdflibs} -o modis_sdev

modis_max : ${OBJS_alert_max}
	${CC} ${OBJS_alert_max} ${hdflibs} -o modis_max

modis_bb : ${OBJS_bb}
	${CC} ${OBJS_bb} ${hdflibs} -o modis_bb

modis_varnti : ${OBJS_varnti}
	${CC} ${OBJS_varnti} ${hdflibs} -o modis_varnti

modis_sdev : ${OBJS_alert_sdev}
	${CC} ${OBJS_alert_sdev} ${hdflibs} -o modis_sdev

makebase : ${OBJS_base}
	${CC} ${OBJS_base} ${hdflibs} -o makebase
tseries : ${OBJS_tseries}
	${CC} ${OBJS_tseries} ${hdflibs} -o tseries


modis_alert : ${OBJS_alert}
	${CC} ${OBJS_alert} ${hdflibs} -o modis_alert

modis_process : ${OBJS}
	${CC} ${OBJS_alert} ${hdflibs} -o modis_process

.cpp.o :
	${CC} ${hdfinc} -c $*.cpp
