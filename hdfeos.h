#ifndef HDFEOS_H
#define HDFEOS_H

//#include "support/configCosmos.h"
//#include "support/datalib.h"
//#include "support/stringlib.h"
#include <mfhdf.h>
#include <hdf.h>
#include <HdfEosDef.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string>
using std::string;
#include <vector>
using std::vector;

struct gctp_file
{
    string path;
    int32_t fileid;
    int32_t gridid;
    uint8_t hindex;
    uint8_t vindex;
    string product;
    uint16_t version;
    uint16_t year;
    uint16_t month;
    uint16_t day;
    uint16_t array[2][1200][1200];
};

struct grid_tile
{
    uint8_t hindex;
    uint8_t vindex;
    float lon_min;
    float lon_max;
    float lat_min;
    float lat_max;
};

int32_t find_tile_paths_system(uint8_t hindex, uint8_t vindex, vector<std::string> &paths);
int32_t find_tile_paths(uint8_t hindex, uint8_t vindex, vector <string> &paths);
int32_t load_tile(string path, gctp_file &tile);
bool data_isdir(string path);
bool data_isfile(string path);
vector < string > string_split(string in, string delimeters);

#endif // HDFEOS_H
