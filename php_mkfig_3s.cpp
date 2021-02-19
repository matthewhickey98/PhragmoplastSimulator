#include <stdlib.h>
#include <string>
#include <iostream>
#include <string.h>
#include "Phragmoplast.h"
#include <stdio.h>
#include <list>
#include <iterator>
#include <Magick++.h>

const char *VERSION = "php_mkfig: v1.0; 10/12/10; B. Piette";
struct Pars
{
    double blur_radius; //  Blur radius
    double grf_max;     //  grf max
    double t_FRAP;      //  grf max
    int no_time;        //  do not write time
} par;

//Parameter Pars[] = {
// {"blur_radius",T_REAL, &par.blur_radius,"0.0",0,"br"," Gaussian Blur Radius"},
// {"grf_max",T_REAL, &par.grf_max,"0.0",0,"M"," Grf max no of crossing"},
// {"t_FRAP",T_REAL, &par.t_FRAP,"0",0,"t"," FRAP time"},
// {"no_time",T_FLAG, &par.no_time,"0",0,"nt"," Do not write time"},
// { 0,T_NONE,0,0,0,0,0},
//};

const char *help[] =
{
    "php_mkfig_3s [-blur_radius=br R] [-grf_max=M R] [-no_time=nt R]",
    0
};

const int PHP_BUFFSIZE = PHP_PATH_LENGTH + 50;

// convert string to a double
double str2dbl(const char *s)
{
    double d;
    sscanf(s, "%lg", &d);
    return (d);
}
// convert string to a int
int str2int(const char *s)
{
    int i;
    sscanf(s, "%d", &i);
    return (i);
}

// BOOST_PYTHON_MODULE(makeFigures)
// {
//     def("runSim", figureMain);
// }


int main(int argc, char **argv)
{
    int i;
    char Buffer[PHP_BUFFSIZE];
    std::list<std::string> files;
    Magick::InitializeMagick(NULL);
    if (argc <= 1)
    {
        std::cerr << "php_mkfig [-blur_radius=br R] [-grf_max=M R] [-no_time=nt R]\n";
        exit(0);
    }

    try
    {
        for (i = 1; i < argc; i++)
        {
            std::cerr << argv[i] << " Radius " << par.blur_radius << "\n";
            if (!strcmp(argv[i], "-blur_radius"))
            {
                par.blur_radius = str2dbl(argv[++i]);
            }
            else if (!strcmp(argv[i], "-grf_max"))
            {
                par.grf_max = str2dbl(argv[++i]);
            }
            else if (!strcmp(argv[i], "-no_time"))
            {
                par.no_time = str2int(argv[++i]);
            }
            else
            {
                files.push_back(std::string(argv[i]));
            }
        }
    }
    catch (std::string s_err)
    {
        std::cerr << "php_mkfig [-blur_radius=br R] [-t_FRAP=t R] [-grf_max=M R] [-no_time=nt R]\n";
        exit(1);
    }

    std::list<std::string>::iterator it;
    for (it = files.begin(); it != files.end(); ++it)
    {
        std::string fname = *it;

        strncpy(Buffer, (*it).c_str(), PHP_BUFFSIZE - 1);
        std::size_t pos = fname.find(".bma");
        if (pos != std::string::npos)
        {
            sprintf(Buffer, "_R%g.jpg", par.blur_radius);
            fname.replace(pos, 4, Buffer);
            std::cerr << "Creating " << fname << "  from " << *it << "\n";

            Bitmap php_bm((*it).c_str());
            php_bm.set_Max(par.grf_max);
            if (par.blur_radius > 0.2)
            {
                php_bm.gaussian_blur(par.blur_radius);
            }
            php_bm.output(fname.c_str(), !par.no_time);
        }
        else
        {
            std::cerr << "Ignoring file " << fname << ": not a .bma file\n";
        }
        ++i;
    }
    std::cerr << "Radius " << par.blur_radius << "\n";
}
