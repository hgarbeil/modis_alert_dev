#include "hdfeos.h"
//#include "support/elapsedtime.h"

vector < vector < grid_tile > > grid_coordinates;

int32_t find_tile_paths_system(uint8_t hindex, uint8_t vindex, vector<std::string> &paths)
{
    char command[200];

    sprintf(command, "/bin/find . -name \"*h%02hu*v%02hu*hdf\"", hindex, vindex);
    FILE *fp = popen(command, "r");

    while ((fgets(command, 200, fp)) != NULL)
    {
        command[strlen(command)-1] = 0;
        paths.push_back(command);
    }
    pclose(fp);
    return paths.size();

}

int32_t find_tile_paths(uint8_t hindex, uint8_t vindex, vector <string> &paths)
{
    struct dir_level
    {
        DIR * jdp;
        dirent * td;
        string path;
        string name;
        double totaltime;
    };

    vector < dir_level > dlevel;
    dir_level tlevel;
    tlevel.path = ".";
    string cpath = ".";
    string cname;
    char gname[6];

    sprintf(gname, "h%02huv%02hu", hindex, vindex);
    paths.clear();

    if ((tlevel.jdp = opendir(tlevel.path.c_str())) != nullptr)
    {
        dlevel.push_back(tlevel);
        do
        {
            while ((dlevel[dlevel.size()-1].td = readdir(dlevel[dlevel.size()-1].jdp)) != nullptr)
            {
//                ElapsedTime et;
                cname = dlevel[dlevel.size()-1].td->d_name;
                cpath = dlevel[dlevel.size()-1].path + "/" + cname;
                if (cname.size() == 45 && cname.find(gname, 17) != string::npos && cname.find(".hdf", 41) != string::npos)
                {
                    paths.push_back(cpath);
                }
                else
                {
                    if (data_isdir(cpath.c_str()))
                    {
                        if (cname[0] != '.')
                        {
                            if ((tlevel.jdp = opendir(cpath.c_str())) != nullptr)
                            {
                                tlevel.path = cpath;
//                                dlevel[dlevel.size()-1].totaltime += et.lap();
                                dlevel.push_back(tlevel);
                                continue;
                            }
                        }
                    }
                }
//                dlevel[dlevel.size()-1].totaltime += et.lap();
            }

            printf("h%u v%u %u %g %u %s\n", hindex, vindex, dlevel.size(), dlevel[dlevel.size()-1].totaltime, paths.size(), dlevel[dlevel.size()-1].path.c_str());
            closedir(dlevel[dlevel.size()-1].jdp);
            dlevel.resize(dlevel.size()-1);
        } while (dlevel.size());
    }
    return paths.size();
}

int32_t load_tile(string path, gctp_file &tile)
{
    int32_t iretn;

    if (data_isfile(path))
    {
        vector < string > path_parts = string_split(path, "/");
        vector < string > name_parts = string_split(path_parts[path_parts.size()-1], ".");

        if (name_parts.size() == 6 && name_parts[5] == "hdf")
        {
            sscanf(path_parts[path_parts.size()-2].data(), "%u.%u.%u", &tile.year, &tile.month, &tile.day);
            vector <string> product_parts = string_split(path_parts[path_parts.size()-3], ".");
            tile.product = product_parts[0];
            tile.version = atoi(product_parts[1].data());
            sscanf(name_parts[2].data(), "h%uv%u", &tile.hindex, &tile.vindex);

            char gridlist[500];
            int32_t ngrid;
            iretn = GDinqgrid((char *)path.c_str(), gridlist, &ngrid);
            if (iretn < 0)
            {
                return iretn;
            }

            tile.fileid = GDopen((char *)path.c_str(), DFACC_READ);
            if (tile.fileid < 0)
            {
                return tile.fileid;
            }

            tile.gridid = GDattach(tile.fileid, gridlist);
            if (tile.gridid < 0)
            {
                return tile.gridid;
            }

            int32 edge[2] = {1200, 1200};

            iretn = GDreadfield(tile.gridid, "LST_Night_1km", NULL, NULL, edge, tile.array[0]);
            if (iretn < 0)
            {
                return iretn;
            }

            iretn = GDreadfield(tile.gridid, "LST_Day_1km", NULL, NULL, edge, tile.array[1]);
            if (iretn < 0)
            {
                return iretn;
            }

            iretn = GDdetach(tile.gridid);
            if (iretn < 0)
            {
                return iretn;
            }

            iretn = GDclose(tile.fileid);
            if (iretn < 0)
            {
                return iretn;
            }

        }

    }
    return 0;
}

bool data_isdir(string path)
{
    struct stat st;

    if (!stat(path.c_str(), &st) && S_ISDIR(st.st_mode))
    {
        return true;
    }
    else
    {
        return false;
    }

}

bool data_isfile(string path)
{
    struct stat st;

    if (!stat(path.c_str(), &st) && S_ISREG(st.st_mode))
    {
        return true;
    }
    else
    {
        return false;
    }

}

vector < string > string_split(string in, string delimeters)
{
    vector<string> result;
    const char *str = in.data();

    do
    {
        const char *begin = str;

        while(*str)
        {
            bool match = false;

            for (size_t i=0; i<delimeters.size(); ++i)
            {
                if (*str == delimeters[i])
                {
                    match = true;
                    break;
                }
            }
            if (match)
            {
                break;
            }
            str++;
        }
        result.push_back(string(begin, str));
    } while (0 != *str++);

    return result;
}
