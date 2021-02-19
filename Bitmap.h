#ifndef Bitmap_H
#define Bitmap_H
#include "php_const.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <Magick++.h>
#include <cmath>
#include <string>
#include <sstream>
using namespace std;
using namespace Magick;
class Bitmap {
public:
    Bitmap(int nx, int ny, double Lx, double Ly, double t, double Max);
    Bitmap(const char *fname);
    ~Bitmap()
    {
        delete array_;
    };

    void clear()
    {
        for (int i = 0; i < (nx_ * ny_); ++i)
        {
            array_[i] = 0.0;
        }
    };

    void line(double x1, double y1, double x2, double y2);

    int x_index_of(double x); // array index of x, -1 -> out of array
    int y_index_of(double y); // array index of y, -1 -> out of array
    inline double x_of(int i)
    {
        return (i * dx_);
    };
    inline double y_of(int j)
    {
        return (j * dy_);
    };
    double &operator()(int ix, int iy)
    {
        return (array_[nx_ * iy + ix]);
    }
    void average(double *v, double Amp); // for kymograph
    double luminosity(double xmin, double ymin,
                      double xmax, double ymax);
    // Compute average luminosity
    void add_line(double *v, int line);                      // add line for kymograph
    double *scan_line(int i)
    {
        return (array_ + i * nx_);
    }; // horizontal line
    void gaussian_blur(double R);

    void output(const char *fname, bool time = false, double t0 = 0); // output graph in a file
    void save(const char *fname);

    inline int nx()
    {
        return (nx_);
    }
    inline int ny()
    {
        return (ny_);
    }
    inline double Lx()
    {
        return (Lx_);
    }
    inline double Ly()
    {
        return (Ly_);
    }
    inline double dx()
    {
        return (dx_);
    }
    inline double dy()
    {
        return (dy_);
    }
    inline void set_t(double t)
    {
        t_ = t;
    }
    inline void set_Max(double Max)
    {
        Max_ = Max;
    }

protected:
    int nx_;
    int ny_;
    double Lx_;
    double Ly_;
    double dx_;
    double dy_;
    double t_;
    double Max_; // Values scalled to G=1
    double *array_;
};

#endif
