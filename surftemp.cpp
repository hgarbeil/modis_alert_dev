#include "surftemp.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
//#include <sys/uio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>

lst_handle lhandle;
gst_handle ghandle;
m02ssh_handle mhandle;



int32_t gst_open()
{
    if ((ghandle.mf = fopen("/local/worldbase/gst/mean_gst","r")) == nullptr)
    {
        return -errno;
    }
    if ((ghandle.sf = fopen("/local/worldbase/gst/std_gst","r")) == nullptr)
    {
        return -errno;
    }
    return 0;
}

int32_t gst_read(vector<lst_coordvalue> &coords)
{
    int32_t iretn;
    for (size_t ic=0; ic<coords.size(); ++ic)
    {
        alert_entry entry;
        entry.latitude = coords[ic].lat * 180. / M_PIl;
        entry.longitude = coords[ic].lon * 180. / M_PIl;
        entry.date = (int)(currentmjd()) + (coords[ic].month - .5) / 12;
        switch (coords[ic].period)
        {
        case modis_period::AQUA_DAY:
            entry.sat = ALERT_SAT_AQUA;
            entry.sun_zenith = 45.;
            break;
        case modis_period::AQUA_NIGHT:
            entry.sat = ALERT_SAT_AQUA;
            entry.sun_zenith = 135.;
            break;
        case modis_period::TERRA_DAY:
            entry.sat = ALERT_SAT_TERRA;
            entry.sun_zenith = 45.;
            break;
        case modis_period::TERRA_NIGHT:
            entry.sat = ALERT_SAT_TERRA;
            entry.sun_zenith = 135.;
            break;
        }
        double result[2];
        iretn = gst_read(entry, result);
        coords[ic].mean = result[0];
        coords[ic].std = result[1];
    }
    return iretn;
}

int32_t gst_read(alert_entry alert, double result[2])
{
    int32_t iretn;
    static uint16_t sat = 5;
    static float sunz = 400.;
    static double date = 0.;
    static double latitude = 100.;
    static double longitude = 100.;
    static size_t soffset = 2488320000L;
    static size_t doffset = 2488320000L;
    static size_t loffset = 2488320000L;
    static size_t offset = 2488320000L;
    static double lresult[2];

    if (ghandle.mf == nullptr || ghandle.sf == nullptr)
    {
        return -200;
    }

    if (alert.sun_zenith != sunz || alert.sat != sat || alert.date != date || alert.longitude != longitude || alert.latitude != latitude)
    {
        if (alert.date != date)
        {
            calstruc cdate = mjd2cal(alert.date);
            doffset = 25920000L * 4 * (cdate.month-1);
            date = alert.date;
        }

        if (alert.sat != sat || alert.sun_zenith != sunz)
        {
            int32_t ihour;
            switch (alert.sat)
            {
            case ALERT_SAT_TERRA:
                // Terra
                if (alert.sun_zenith < 90.)
                {
                    // Day
                    ihour = 1;
                }
                else
                {
                    // Night
                    ihour = 3;
                }
                break;
            case ALERT_SAT_AQUA:
                // Aqua
                if (alert.sun_zenith < 90.)
                {
                    // Day
                    ihour = 2;
                }
                else
                {
                    // Night
                    ihour = 0;
                }
                break;
            }
            soffset = 25920000L * ihour;
            sat = alert.sat;
            sunz = alert.sun_zenith;
        }

        if (alert.longitude != longitude || alert.latitude != latitude)
        {
            loffset = (int32_t)((alert.longitude + 180.025) / .05) + 7200 * (int32_t)((90.025 - alert.latitude) / .05);
            latitude = alert.latitude;
            longitude = alert.longitude;
        }

        offset = 2L * (doffset + soffset + loffset);

        uint16_t input;
        iretn = fseek(ghandle.mf, offset, SEEK_SET);
        if (iretn >= 0)
        {
            iretn = fread(&input, 2, 1, ghandle.mf);
            if (iretn == 1)
            {
                lresult[0] = input / 100.;
            }
        }
        iretn = fseek(ghandle.sf, offset, SEEK_SET);
        if (iretn >= 0)
        {
            iretn = fread(&input, 2, 1, ghandle.sf);
            if (iretn == 1)
            {
                lresult[1] = input / 100.;
            }
        }
    }

    result[0] = lresult[0];
    result[1] = lresult[1];

    return 0;

}

int32_t gst_close()
{
    if (ghandle.mf != nullptr)
    {
        fclose(ghandle.mf);
    }
    if (ghandle.sf != nullptr)
    {
        fclose(ghandle.sf);
    }
    return 0;
}

int32_t lst_open(string base)
{
    lhandle.base = base;
    lhandle.qindex = lst_max_quad;
    lhandle.qcount = 0;
    lhandle.quads.resize(lst_max_quad);
    for (size_t i=0; i<lst_max_quad; ++i)
    {
        lhandle.quads[i].ff = nullptr;
        lhandle.quads[i].mean.resize(lst_size);
        for (size_t j=0; j<lst_size; ++j)
        {
            lhandle.quads[i].mean[j].resize(lst_size);
        }
        lhandle.quads[i].std.resize(lst_size);
        for (size_t j=0; j<lst_size; ++j)
        {
            lhandle.quads[i].std[j].resize(lst_size);
        }
    }
    return 0;
}

int32_t lst_read(vector<lst_coordvalue> &coords)
{
    for (size_t ic=0; ic<coords.size(); ++ic)
    {
        lst_index(coords[ic]);
        if (lhandle.qcount == 0 || coords[ic].vgrid != lhandle.vgrid || coords[ic].hgrid != lhandle.hgrid)
        {
            lst_load(coords[ic]);
        }
        coords[ic].mean = lhandle.quads[lhandle.qindex].mean[coords[ic].vindex][coords[ic].hindex];
        coords[ic].std = lhandle.quads[lhandle.qindex].std[coords[ic].vindex][coords[ic].hindex];
    }

    return 0;

}

int32_t lst_load(lst_coordvalue &coord)
{
    if (lhandle.qcount && lhandle.hgrid == coord.hgrid && lhandle.vgrid == coord.vgrid)
    {
        return 0;
    }

    for (size_t iq=0; iq<lhandle.qcount; ++iq)
    {
        if (lhandle.quads[iq].hgrid == coord.hgrid && lhandle.quads[iq].vgrid == coord.vgrid)
        {
            lhandle.hgrid = lhandle.quads[iq].hgrid;
            lhandle.vgrid = lhandle.quads[iq].vgrid;
            lhandle.qindex = iq;
            return 0;
        }
    }

    if (lhandle.qcount < lst_max_quad)
    {
        lhandle.qindex = lhandle.qcount++;
    }
    else
    {

        double oldestt = currentmjd();
        size_t oldesti = lst_max_quad+1;
        for (size_t i=0; i<lhandle.qcount; ++i)
        {
            if (lhandle.quads[i].date < oldestt)
            {
                oldestt = lhandle.quads[i].date;
                oldesti = i;
            }
        }
        if (lhandle.quads[oldesti].ff != nullptr)
        {
            fclose(lhandle.quads[oldesti].ff);
        }
        lhandle.qindex = oldesti;
    }
    lhandle.hgrid = coord.hgrid;
    lhandle.vgrid = coord.vgrid;

    char path[200];
    sprintf(path, "/local/worldbase/%s/%02u_%1u/h%02uv%02u.bsq", lhandle.base.c_str(), coord.month, coord.period, coord.hgrid, coord.vgrid);
    FILE *tff = fopen(path, "rb");
    if (tff == nullptr)
    {
        return -errno;
    }
    lhandle.quads[lhandle.qindex].ff = tff;

    lhandle.quads[lhandle.qindex].month = coord.month;
    lhandle.quads[lhandle.qindex].period = coord.period;
    lhandle.quads[lhandle.qindex].vgrid = coord.vgrid;
    lhandle.quads[lhandle.qindex].hgrid = coord.hgrid;

    for (size_t i=0; i<lst_size; ++i)
    {
        if (lhandle.quads[lhandle.qindex].ff == nullptr)
        {
            memset(lhandle.quads[lhandle.qindex].mean[i].data(), 0, sizeof(float) * lst_size);
        }
        else
        {
            fread(lhandle.quads[lhandle.qindex].mean[i].data(), sizeof(float), lst_size, lhandle.quads[lhandle.qindex].ff);
        }
    }

    for (size_t i=0; i<lst_size; ++i)
    {
        if (lhandle.quads[lhandle.qindex].ff == nullptr)
        {
            memset(lhandle.quads[lhandle.qindex].std[i].data(), 0, sizeof(float) * lst_size);
        }
        else
        {
            fread(lhandle.quads[lhandle.qindex].std[i].data(), sizeof(float), lst_size, lhandle.quads[lhandle.qindex].ff);
        }
    }
    lhandle.quads[lhandle.qindex].date = currentmjd();


    return 0;
}

int32_t lst_close()
{
    for (size_t iq=0; iq<lhandle.qcount; ++iq)
    {
        fclose(lhandle.quads[iq].ff);
    }
    return 0;
}

int32_t lst_index(lst_coordvalue &coord)
{
    double lat = fmod(coord.lat + M_PI_2l, M_PIl) - M_PI_2l;
    double slat = lat / lst_step;
    coord.vgrid = 8 - floor(slat);
    if (coord.vgrid > 17)
    {
        coord.vgrid = 0;
    }
    coord.vindex = -(lst_size * (slat - (9 - (int16_t)coord.vgrid))) - .5;

    double lon = fmod(coord.lon + M_PIl, 2. * M_PIl) - M_PIl;
    double slon = cos(lat) * lon / lst_step;
    coord.hgrid = 18 + floor(slon);
    if (coord.hgrid > 35)
    {
        coord.hgrid = 0;
    }
    coord.hindex = lst_size * (slon - ((int16_t)coord.hgrid - 18)) - .5;

    return 0;
}

int32_t lst_coord(lst_coordvalue &coord)
{
    coord.lat = lst_step * ((9 - (int16_t)coord.vgrid) - (coord.vindex + .5) / lst_size);
    coord.lon = (lst_step / cos(coord.lat)) * (((int16_t)coord.hgrid - 18) + (coord.hindex + .5) / lst_size);
    return 0;
}

int32_t m02ssh_index(m02ssh_coordvalue &coord)
{
    double lat = fmod(coord.lat + M_PI_2, M_PI) - M_PI_2;
    coord.vindex = (M_PI_2 - lat) / m02ssh_step;

    double lon = fmod(coord.lon + M_PI, 2. * M_PI) - M_PI;
    coord.hindex = m02ssh_height + .5 + (cos(lat) * lon) / m02ssh_step;

    return 0;
}

int32_t m02ssh_coord(m02ssh_coordvalue &coord)
{
    coord.lat = M_PI_2 - m02ssh_step * (coord.vindex + .5);
    coord.lon = (m02ssh_step / cos(coord.lat)) * (coord.hindex + .5);
    return 0;
}

calstruc mjd2cal(double mjd)
{
    static double lmjd = 0.;
    static calstruc date;

    if (lmjd != mjd)
    {
        double dom;
        double doy;

        lmjd = mjd;

        mjd2ymd(mjd, date.year, date.month, dom, doy);
        date.doy = (int32_t)doy;
        date.dom = (int32_t)dom;
        doy = (doy - date.doy) * 24.;
        date.hour = (int32_t)doy;
        doy = (doy - date.hour) * 60.;
        date.minute = (int32_t)doy;
        doy = (doy - date.minute) * 60.;
        date.second = (int32_t)doy;
        doy = (doy - date.second) * 1e9;
        date.nsecond = (int32_t)(doy + .5);
    }

    return date;
}

int32_t m02ssh_open(string product, uint16_t jday, modis_period period)
{
    char path[200];

    if (product == "count")
    {
        mhandle.count.resize(m02ssh_height);

        sprintf(path, "/local/worldbase/02ssh/count_%03u_%02hhu.bsq", jday, static_cast<uint8_t>(period));
        FILE *cfp = fopen(path, "rb");
        if (cfp == nullptr)
        {
            return -errno;
        }

        for (uint16_t ir=0; ir<m02ssh_height; ++ir)
        {
            mhandle.count[ir].resize(m02ssh_width);
            fread(mhandle.count[ir].data(), sizeof(uint16_t), m02ssh_width, cfp);
        }

        fclose(cfp);
    }
    else
    {
        mhandle.mean.resize(m02ssh_height);
        mhandle.std.resize(m02ssh_height);
        // add v2 to get the mean value - rather than upper quadrile
        //sprintf(path, "/local/worldbase/02ssh/v2/%s_%03u_%02hhu.bsq", product.c_str(), jday, static_cast<uint8_t>(period));
        sprintf(path, "/local/worldbase/02ssh/%s_%03u_%02hhu.bsq", product.c_str(), jday, static_cast<uint8_t>(period));
        FILE *dfp = fopen(path, "rb");
        if (dfp == nullptr)
        {
            return -errno;
        }

        for (uint16_t ir=0; ir<m02ssh_height; ++ir)
        {
            mhandle.mean[ir].resize(m02ssh_width);
            fread(mhandle.mean[ir].data(), sizeof(float), m02ssh_width, dfp);
        }

        for (uint16_t ir=0; ir<m02ssh_height; ++ir)
        {
            mhandle.std[ir].resize(m02ssh_width);
            fread(mhandle.std[ir].data(), sizeof(float), m02ssh_width, dfp);
        }

        fclose(dfp);
    }

    return 0;
}

int32_t m02ssh_read(vector<m02ssh_coordvalue> &coords)
{
    for (size_t ic=0; ic<coords.size(); ++ic)
    {
        m02ssh_index(coords[ic]);
        if (mhandle.mean.size())
        {
            coords[ic].mean = mhandle.mean[coords[ic].vindex][coords[ic].hindex];
        }
        if (mhandle.std.size())
        {
            coords[ic].std = mhandle.std[coords[ic].vindex][coords[ic].hindex];
        }
        if (mhandle.count.size())
        {
            coords[ic].count = mhandle.count[coords[ic].vindex][coords[ic].hindex];
        }
    }

    return 0;

}

int32_t mjd2ymd(double mjd, int32_t &year, int32_t &month, double &day, double &doy)
{
    static double lmjd = 0.;
    static int32_t lyear = 1858;
    static int32_t lmonth = 11;
    static double lday = 17.;
    static double ldoy = 321.;

    if (mjd != lmjd)
    {
        int32_t a, b, c, d, e, z, alpha;
        double f;

        lmjd = mjd;
        mjd += 2400001.;
        z = (int32_t)mjd;
        f = mjd - z;

        if (z<2299161)
            a = z;
        else
        {
            alpha = (int32_t)((z - 1867216.25)/36524.25);
            a = z +1 + alpha - (int32_t)(alpha/4);
        }

        b = a + 1524;
        c = (int32_t)((b - 122.1)/365.25);
        d = (int32_t)(365.25*c);
        e = (int32_t)((b - d)/30.6001);

        lday = b - d - (int32_t)(30.6001 * e) + f;
        if (e < 14)
            lmonth = e - 1;
        else
            lmonth = e - 13;
        if (lmonth > 2)
            lyear = c - 4716;
        else
            lyear = c - 4715;
        ldoy = (int32_t)((275 * lmonth)/9) - (2-isleap(lyear))*(int32_t)((lmonth+9)/12) + lday - 30;
    }

    year = lyear;
    month = lmonth;
    day = lday;
    doy  = ldoy;
    return 0;
}

double currentmjd(double offset)
{
    double mjd;

    // unfortunatelly MSVC does not support gettimeofday
#ifdef COSMOS_WIN_BUILD_MSVC
    TimeUtils tu;
    mjd = unix2utc(tu.secondsSinceEpoch() + _timezone);
#else
    struct timeval mytime;
    gettimeofday(&mytime, NULL);
    mjd = unix2utc(mytime);
#endif
    return mjd+offset;
}

double currentmjd()
{
    return currentmjd(0.);
}

int16_t isleap(int32_t year)
{
    if (!(year % 4))
    {
        if (!(year%100))
        {
            if (!(year%400))
            {
                return (1);
            }
            else
            {
                return 0;
            }
        }
        else
        {
            return 1;
        }
    }
    else
    {
        return 0;
    }
}

double unix2utc(struct timeval unixtime)
{
    double utc;
    utc = 40587. + (unixtime.tv_sec + unixtime.tv_usec / 1000000.) / 86400.;

    return utc;
}

float floatfrom(uint8_t *pointer, ByteOrder order)
{
    float result;
    uint8_t *rb;

    rb = (uint8_t *)&result;
    if (local_byte_order() == order)
    {
        memcpy((void *)rb,pointer,4);
    }
    else
    {
        rb[3] = pointer[0];
        rb[2] = pointer[1];
        rb[1] = pointer[2];
        rb[0] = pointer[3];
    }

    return (result);
}

ByteOrder local_byte_order()
{
    uint16_t test = 1;
    uint8_t *check;

    check = (uint8_t *)&test;

    if (check[0] == 0)
        return (ByteOrder::BIGENDIAN);
    else
        return (ByteOrder::LITTLEENDIAN);
}

