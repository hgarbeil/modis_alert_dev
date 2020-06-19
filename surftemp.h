#ifndef ALERTDB_H
#define ALERTDB_H 1
#endif
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "trig.h"

#define TEMPB21(x)   (3634.22/log(1.2246e11/(1e6*x)+1.))
#define TEMP32(x) (1197.00/log(474.651/(x)+1.)) 
#define ALERT_SAT_TERRA			1
#define ALERT_SAT_AQUA			2

enum class ByteOrder : std::uint8_t {
    //! Big Endian byte order
    BIGENDIAN=0, // was previouly ORDER_BIGENDIAN, replace by ByteOrder::BIGENDIAN
    //! PowerPC byte order
    PPC=ByteOrder::BIGENDIAN,
    //! Motorola byte order
    MOTOROLA=ByteOrder::BIGENDIAN,
    //! Little Endian byte order
    LITTLEENDIAN=1, // was previouly ORDER_LITTLEENDIAN
    //! Intel byte order
    INTEL=ByteOrder::LITTLEENDIAN,
    //! Network byte order
    NETWORK=ByteOrder::BIGENDIAN
};


struct alert_entry
    {
    double date;
    double longitude;
    double latitude;
    uint16_t sat;
    float modis21;
    float modis22;
    float modis6;
    float modis31;
    float modis32;
    float ratio;
    uint16_t line;
    uint16_t sample;
    float zenith;
    float azimuth;
    float sun_zenith;
    float sun_azimuth;
    float glint_angle;
    int32_t next_entry;
    };

struct gst_handle
{
	FILE *mf;
	FILE *sf;
} ;

const size_t lst_size = 1200;
const double modis_radius = 6371007.181;
const double lst_step = M_PI / 18.;
const size_t lst_max_quad = 36;

const uint16_t m02ssh_height = 240 * 18;
const uint16_t m02ssh_width = 2 * m02ssh_height;
const double m02ssh_step = M_PI / (m02ssh_height);

enum class modis_period : size_t
{
    AQUA_NIGHT,
    TERRA_DAY,
    AQUA_DAY,
    TERRA_NIGHT
};

struct m02ssh_coordvalue
{
    size_t vindex;
    size_t hindex;
    float lat;
    float lon;
    uint16_t jday;
    modis_period period;
    float mean;
    float std;
    uint16_t count;
};

struct m02ssh_handle
{
    bool loaded = false;
    vector<vector<float> > mean;
    vector<vector<float> > std;
    vector<vector<uint16_t>> count;
};

struct lst_coordvalue
{
    size_t vgrid;
    size_t hgrid;
    size_t vindex;
    size_t hindex;
    double lat;
    double lon;
    size_t month;
    modis_period period;
    double mean;
    double std;
};

struct lst_quad
{
    size_t vgrid;
    size_t hgrid;
    size_t month;
    modis_period period;
    double date;
    FILE *ff;
    vector<vector<float> > mean;
    vector<vector<float> > std;
} ;

struct lst_handle
{
    string base;
    size_t hgrid;
    size_t vgrid;
    size_t qindex;
    size_t qcount;
    vector<lst_quad> quads;
} ;

struct calstruc
{
    int32_t year;
    int32_t month;
    int32_t dom;
    int32_t doy;
    int32_t hour;
    int32_t minute;
    int32_t second;
    int32_t nsecond;
};


int32_t m02ssh_index(m02ssh_coordvalue &coord);
int32_t mo2ssh_coord(m02ssh_coordvalue &coord);
int32_t m02ssh_open(string product, uint16_t jday, modis_period period);
int32_t m02ssh_read(vector<m02ssh_coordvalue> &coords);

int32_t gst_read(vector<lst_coordvalue> &coords);
int32_t gst_read(alert_entry alert, double result[]);
int32_t gst_open();
int32_t gst_close();

int32_t lst_read(vector<lst_coordvalue> &coords);
int32_t lst_open(string base="lst");
int32_t lst_load(lst_coordvalue &coord);
int32_t lst_close();
int32_t lst_index(lst_coordvalue &coord);
int32_t lst_coord(lst_coordvalue &coord);
calstruc mjd2cal(double mjd);
int32_t mjd2ymd(double mjd, int32_t &year, int32_t &month, double &day, double &doy);
double currentmjd(double offset);
double currentmjd();
int16_t isleap(int32_t year);
double unix2utc(struct timeval unixtime);
float floatfrom(uint8_t *pointer, ByteOrder order);
ByteOrder local_byte_order();

//#ifdef __cplusplus
//}
//#endif

