#ifndef Bitmap_HPP
#define Bitmap_HPP

#include "php_const.hpp"

using namespace std;
// using namespace Magick;

class Bitmap
{
protected:
    int numX;
    int numY;
    double lenX;
    double lenY;
    double deltaX;
    double deltaY;
    double time;
    double maxValue; // Values scalled to G=1
    std::vector<double> imageArray;

public:
    Bitmap()
    {
    }
    Bitmap(int nx, int ny, double Lx, double Ly, double t, double Max);
    Bitmap(const std::string &fname);
    ~Bitmap(){};
    // void SaveAsMatFile( const std::string& name) const;
    void Clear()
    {
        imageArray.clear();
        imageArray.resize(numX * numY, 0);
    };

    void Line(double x1, double y1, double x2, double y2);

    int IndexOfX(double x); // array index of x, -1 -> out of array
    int IndexOfY(double y); // array index of y, -1 -> out of array
    inline double x_of(int i)
    {
        return (i * deltaX);
    };
    inline double y_of(int j)
    {
        return (j * deltaY);
    };
    double &operator()(int ix, int iy)
    {
        return (imageArray[numX * iy + ix]);
    }
    void Average(std::vector<double> &v, double Amp);                      // for kymograph
    double Luminosity(double xmin, double ymin, double xmax, double ymax); // Compute average luminosity

    void AddLine(std::vector<double> &v, int line); // add line for kymograph

    double &Scan_line(int i)
    {
        return (imageArray[i * numX]);
    }; // horizontal line

    void GaussianBlur(double R);

    void Output(const std::string &fname, bool writeTime, double t0); 
    void Save(const std::string &fname);

    inline int nx() const
    {
        return (numX);
    }
    inline int ny() const
    {
        return (numY);
    }
    inline double Lx()
    {
        return (lenX);
    }
    inline double Ly()
    {
        return (lenY);
    }
    inline double dx()
    {
        return (deltaX);
    }
    inline double dy()
    {
        return (deltaY);
    }
    inline void set_t(double t)
    {
        time = t;
    }
    inline void set_Max(double Max)
    {
        maxValue = Max;
    }
};

#endif
