#include "Bitmap.hpp"
#include "Magick++.h"
#include <functional>

/**
 * @brief Construct a new Bitmap for use within the simulation
 * 
 * @param nx    number of x pixels
 * @param ny    number of y pixels
 * @param Lx    length of x axis
 * @param Ly    length of y axis
 * @param t     time
 * @param Max   maximum superposition
 */
Bitmap::Bitmap(int nx, int ny, double Lx, double Ly, double t, double Max)
    : numX(nx), numY(ny), lenX(Lx), lenY(Ly), time(t), maxValue(Max)
{
    imageArray.resize(nx * ny, 0);
    deltaX = lenX / (numX - 1.0);
    deltaY = lenY / (numY - 1.0);
    Clear();
}

/**
 * @brief Construct a new Bitmap from a file
 * 
 * @param fname     name of file
 */
Bitmap::Bitmap(const std::string &fname) : numX(0), numY(0)
{
    ifstream ifs(fname, ofstream::binary);
    ifs.seekg(0);

    ifs.read((char *)&numX, sizeof(numX));
    ifs.read((char *)&numY, sizeof(numY));
    ifs.read((char *)&lenX, sizeof(lenX));
    ifs.read((char *)&lenY, sizeof(lenY));
    ifs.read((char *)&time, sizeof(time));
    deltaX = lenX / (numX - 1);
    deltaY = lenY / (numY - 1);

    imageArray.resize(numX * numY, 0);
    ifs.read((char *)(&imageArray[0]), sizeof(double) * numX * numY);
}

/**
 * @brief Add a line to the bitmap
 * 
 * @param x1    x coordinate of first point
 * @param y1    y coordinate of first point
 * @param x2    x coordinate of second point
 * @param y2    y coordinate of second point
 */
void Bitmap::Line(double x1, double y1, double x2, double y2)
{
    int i1, i2, i, j;
    double x, y;

    i1 = IndexOfX(x1);
    i2 = IndexOfX(x2);
    if (i1 < 0)
    {
        i1 = (x1 < 0) ? 0 : numX - 1;
    } //{ i1 = 0; }
    if (i2 < 0)
    {
        i2 = (x2 < 0) ? 0 : numX - 1;
    } //{ i2 = 0; }
    if (i1 >= numX)
    {
        i1 = numX - 1;
    }
    if (i2 >= numX)
    {
        i2 = numX - 1;
    }

    if (i2 < i1)
    {
        i = i2;
        i2 = i1;
        i1 = i;
    }

    // cout<<"i1="<<i1<<"  i2="<<i2<<"\n";
    if ((i1 >= 0) && (i2 >= 0))
    {
        for (i = i1; i <= i2; ++i)
        {
            x = x_of(i);
            y = ((y2 - y1) * x - y2 * x1 + y1 * x2) / (x2 - x1);
            j = IndexOfY(y);
            if (j >= 0)
            {
                (*this)(i, j) += 1;
            }
        }
    }
}

/**
 * @brief Get the index of the x coordinate
 * 
 * @param x     x coordinate
 * @return int  index
 */
int Bitmap::IndexOfX(double x)
{
    if ((x < 0) || (x > lenX))
    {
        return (-1);
    } // out of bound
    return ((int)floor(x / lenX * (numX - 1.0) + 0.5));
}

/**
 * @brief Get the index of the y coordinate
 * 
 * @param y     y coordinate
 * @return int  index
 */
int Bitmap::IndexOfY(double y)
{
    if ((y < 0) || (y > lenY))
    {
        return (-1);
    } // out of bound
    return ((int)floor(y / lenY * (numY - 1.0) + 0.5));
}

/**
 * @brief Compute the vertical average value for a kymograph
 * 
 * @param v     vector to store results
 * @param Amp   amplitude
 */
void Bitmap::Average(std::vector<double> &v, double Amp)
{
    int i, j;

    if (Amp <= 0)
    {
        Amp = 1;
    }

    for (i = 0; i < numX; ++i)
    {
        v[i] = 0;
        for (j = 0; j < numY; ++j)
        {
            v[i] += imageArray[i + numX * j];
        }
        v[i] *= Amp / numY;
    }
}


/**
 * @brief Compute the average total luminosity inside the specified region
 * 
 * @param xmin      x coordinate of the lower left corner
 * @param ymin      y coordinate of the lower left corner
 * @param xmax      x coordinate of the upper right corner
 * @param ymax      y coordinate of the upper right corner
 * @return double   average total luminosity
 */
double Bitmap::Luminosity(double xmin, double ymin, double xmax, double ymax)
{
    int i, j, n;
    double x, y, v;

    v = 0;
    n = 0;
    cout << "t=" << time << "\n";
    for (i = 0; i < numX; ++i)
    {
        x = x_of(i);
        if ((xmin <= x) && (x <= xmax))
        {
            for (j = 0; j < numY; ++j)
            {
                y = y_of(j);
                if ((ymin <= y) && (y <= ymax))
                {
                    v += imageArray[i + numX * j];
                    // if(array_[i + numX*j]) { cout<<"a("<<i<<","<<j<<")="<<array_[i + numX*j]<<"\n";}
                    ++n;
                }
            }
        }
    }
    return (v);
}

#ifdef OLDSTFF
// Compute the average total luminosity inside the specified region
double Bitmap::luminosity(double xmin, double ymin, double xmax, double ymax)
{
    int i, j, n, imin, imax, jmin, jmax;
    double v;

    v = 0;
    n = 0;
    // Make sure we take the inside of the FRAp region and in a symetric way
    imin = IndexOfX(xmin);
    imax = IndexOfX(xmax);
    jmin = IndexOfY(ymin);
    jmax = IndexOfY(ymax);

    for (i = imin; i <= imax; ++i)
    {
        for (j = jmin; j <= jmax; ++j)
        {
            v += array_[i + numX * j];
            ++n;
        }
    }
    return (v);
}
#endif

/**
 * @brief Add a line to the kymograph
 * 
 * @param v     vector to store results
 * @param line  line to add
 */
void Bitmap::AddLine(std::vector<double> &v, int line)
{
    if (line < numY)
    {
        for (int i = 0; i < numX; ++i)
        {
            imageArray[numX * line + i] = v[i];
        }
    }
}

/**
 * @brief Perform a gaussian blur on the bitmap
 * 
 * @param R     radius of the gaussian blur
 */
void Bitmap::GaussianBlur(double R)
{
    int i, j, li, lj, ni, i0, j0;
    double r;

    double norm;
    R *= R;
    norm = 1.0 / (M_PI * R);
    // cout<<"gaussian_blur()\n";
    std::vector<double> larray(numX * numY);

    if (R > 0.2)
    {
        ni = (int)floor(R * 4);
        for (j0 = 0; j0 < numY; ++j0)
        {
            for (i0 = 0; i0 < numX; ++i0)
            {
                larray[i0 + j0 * numX] = 0;
                for (lj = -ni; lj <= ni; ++lj)
                {
                    j = j0 + lj;
                    if ((j >= 0) && (j < numY))
                    {
                        for (li = -ni; li <= ni; ++li)
                        {
                            i = i0 + li;
                            if ((i >= 0) && (i < numX))
                            {
                                r = sqrt(li * li + lj * lj);
                                if (r <= R)
                                {
                                    larray[i0 + j0 * numX] += imageArray[i + j * numX] * exp(-r * r / R) * norm;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    imageArray = larray;
}

// Output figure into file
// time : if true write time on fig
/**
 * @brief Save the bitmap figure to a JPEG file
 * 
 * @param fname         name of the file
 * @param writeTime     write time on the figure
 * @param t0            time of the first frame
 */
void Bitmap::Output(const std::string &fname, bool writeTime, double t0)
{
    char dimbuff[PHP_PATH_LENGTH];
    double RGB_array[numX * numY * 3];
    double v, vmax = 0;
    int i, j;
    Magick::ColorRGB c_red(1, 0, 0);

    if (maxValue <= 0)
    {
        vmax = 0;
        i = j = 0;
        v = *(std::max_element(imageArray.begin(), imageArray.end()));
        if (v > vmax)
        {
            vmax = v;
        }
        cout << "Vmax from array is: " << vmax << endl;
    }
    else
    {
        vmax = maxValue;
    }

    i = j = 0;
    cout << "max(imageArray) = " << *(std::max_element(imageArray.begin(), imageArray.end())) << endl;
    cout << "NumX: " << numX << "\tNumY: " << numY << endl;

    while (i < numX * numY)
    {
        RGB_array[j++] = 0;             // R
        RGB_array[j++] = imageArray[i]; // G
        RGB_array[j++] = 0;             // B
        i++;
    }
    if (maxValue <= 0)
        std::cout << "vmax=" << vmax << "\n";

    stringstream ss;
    ss << "Image_T_" << time;

    sprintf(dimbuff, "%dx%d", numX, numY);

    Magick::Image *image;

    image = new Magick::Image(dimbuff, "white");
    image->read(numX, numY, "RGB", Magick::DoublePixel, RGB_array);

    // Draw time on figure

    if (writeTime)
    {
        const string font = "Helvetica";

        stringstream ss;
        ss << (time - t0);
        string t_buff = ss.str();

        Magick::DrawableText dt(6, 22, t_buff);

        image->fontPointsize(16);
        image->font(font);
        image->strokeColor(c_red);
        image->strokeWidth(1);
        image->draw(dt);
    }

    string newFname(fname);
    image->write(newFname);
    // cout << "Wrote image: " << newFname << endl;
    cout << "time: " << time << "\tt: " << t0 << endl;

    // delete image;
}

/**
 * @brief Save the bitmap figure to a BMP file
 * 
 * @param fname     name of the file
 */
void Bitmap::Save(const std::string &fname)
{
    ofstream ofs(fname, ofstream::binary);

    ofs.write((char *)&numX, sizeof(numX));
    ofs.write((char *)&numY, sizeof(numY));
    ofs.write((char *)&lenX, sizeof(lenX));
    ofs.write((char *)&lenY, sizeof(lenY));
    ofs.write((char *)&time, sizeof(time));
    ofs.write((char *)(&imageArray[0]), sizeof(double) * numX * numY);

    ofs.close();
}
