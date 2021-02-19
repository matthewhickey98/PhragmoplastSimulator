#include "Bitmap.h"

// class php_bitmap
Bitmap::Bitmap(int nx, int ny, double Lx, double Ly, double t,
               double Max) : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly), t_(t), Max_(Max)
{
    array_ = new double[nx_ * ny_];
    for (int i = 0; i < nx_ * ny_; ++i)
    {
        array_[i] = 0.0;
    }
    dx_ = Lx_ / (nx_ - 1.0);
    dy_ = Ly_ / (ny_ - 1.0);

    clear();
}

// Create php_bitmap from file data
Bitmap::Bitmap(const char *fname) : nx_(0), ny_(0)
{
    ifstream ifs(fname, ofstream::binary);
    ifs.seekg(0);

    ifs.read((char *)&nx_, sizeof(nx_));
    ifs.read((char *)&ny_, sizeof(ny_));
    ifs.read((char *)&Lx_, sizeof(Lx_));
    ifs.read((char *)&Ly_, sizeof(Ly_));
    ifs.read((char *)&t_, sizeof(t_));

    dx_ = Lx_ / (nx_ - 1);
    dy_ = Ly_ / (ny_ - 1);

    //  std::cerr<<"nx_="<<nx_<<" ny_="<<ny_<<" Lx_="<<Lx_<<" Ly_="<<Ly_<<"\n";
    array_ = new double[nx_ * ny_];

    ifs.read((char *)array_, sizeof(double) * nx_ * ny_);
}

void Bitmap::line(double x1, double y1, double x2, double y2)
{
    int i1, i2, i, j;
    double x, y;

    i1 = x_index_of(x1);
    i2 = x_index_of(x2);
    if (i1 < 0)
    {
        i1 = (x1 < 0) ? 0 : nx_ - 1;
    } //{ i1 = 0; }
    if (i2 < 0)
    {
        i2 = (x2 < 0) ? 0 : nx_ - 1;
    } //{ i2 = 0; }
    if (i1 >= nx_)
    {
        i1 = nx_ - 1;
    }
    if (i2 >= nx_)
    {
        i2 = nx_ - 1;
    }

    if (i2 < i1)
    {
        i = i2;
        i2 = i1;
        i1 = i;
    }

    //cerr<<"i1="<<i1<<"  i2="<<i2<<"\n";
    if ((i1 >= 0) && (i2 >= 0))
    {
        for (i = i1; i <= i2; ++i)
        {
            x = x_of(i);
            y = ((y2 - y1) * x - y2 * x1 + y1 * x2) / (x2 - x1);
            j = y_index_of(y);
            if (j >= 0)
            {
                (*this)(i, j) += 1;
            }
            //cerr<<"i="<<i<<" j="<<j<<" x="<<x<<" y="<<y<<"\n";
        }
    }
}

int Bitmap::x_index_of(double x)
{
    if ((x < 0) || (x > Lx_))
    {
        return (-1);
    } // out of bound
    return ((int)floor(x / Lx_ * (nx_ - 1.0) + 0.5));
}

int Bitmap::y_index_of(double y)
{
    if ((y < 0) || (y > Ly_))
    {
        return (-1);
    } // out of bound
    return ((int)floor(y / Ly_ * (ny_ - 1.0) + 0.5));
}

// Compute verical average value for kymogrpah
// return results in v
void Bitmap::average(double *v, double Amp)
{
    int i, j;

    if (Amp <= 0)
    {
        Amp = 1;
    }

    for (i = 0; i < nx_; ++i)
    {
        v[i] = 0;
        for (j = 0; j < ny_; ++j)
        {
            v[i] += array_[i + nx_ * j];
        }
        v[i] *= Amp / ny_;
    }
}

// Compute the average total luminosity inside the specified region
double Bitmap::luminosity(double xmin, double ymin,
                          double xmax, double ymax)
{
    int i, j, n;
    double x, y, v;

    v = 0;
    n = 0;
    cerr << "t=" << t_ << "\n";
    for (i = 0; i < nx_; ++i)
    {
        x = x_of(i);
        if ((xmin <= x) && (x <= xmax))
        {
            for (j = 0; j < ny_; ++j)
            {
                y = y_of(j);
                if ((ymin <= y) && (y <= ymax))
                {
                    v += array_[i + nx_ * j];
                    //if(array_[i + nx_*j]) { cerr<<"a("<<i<<","<<j<<")="<<array_[i + nx_*j]<<"\n";}
                    ++n;
                }
            }
        }
    }
    //cerr<<"n="<<n<<"\n";
    return (v);
}

#ifdef OLDSTFF
// Compute the average total luminosity inside the specified region
double Bitmap::luminosity(double xmin, double ymin,
                          double xmax, double ymax)
{
    int i, j, n, imin, imax, jmin, jmax;
    double v;

    v = 0;
    n = 0;
    // Make sure we take the inside of the FRAp region and in a symetric way
    imin = x_index_of(xmin);
    imax = x_index_of(xmax);
    jmin = y_index_of(ymin);
    jmax = y_index_of(ymax);

    for (i = imin; i <= imax; ++i)
    {
        for (j = jmin; j <= jmax; ++j)
        {
            v += array_[i + nx_ * j];
            ++n;
        }
    }
    return (v);
}
#endif

// Add one line to kymograph
void Bitmap::add_line(double *v, int line)
{
    if (line < ny_)
    {
        for (int i = 0; i < nx_; ++i)
        {
            array_[nx_ * line + i] = v[i];
        }
    }
}

// Perform a gaussian blur with radius R
// int_plane exp(-r^2/R^2) rdr dphi = pi R^2
void Bitmap::gaussian_blur(double R)
{
    int i, j, li, lj, ni, i0, j0;
    double r;
    double *larray, norm;
    R *= R;
    norm = 1.0 / (M_PI * R);
    //cerr<<"gaussian_blur()\n";
    larray = new double[nx_ * ny_];

    if (R > 0.2)
    {
        ni = (int)floor(R * 4);
        for (j0 = 0; j0 < ny_; ++j0)
        {
            for (i0 = 0; i0 < nx_; ++i0)
            {
                larray[i0 + j0 * nx_] = 0;
                for (lj = -ni; lj <= ni; ++lj)
                {
                    j = j0 + lj;
                    if ((j >= 0) && (j < ny_))
                    {
                        for (li = -ni; li <= ni; ++li)
                        {
                            i = i0 + li;
                            if ((i >= 0) && (i < nx_))
                            {
                                r = sqrt(li * li + lj * lj);
                                if (r <= R)
                                {
                                    larray[i0 + j0 * nx_] += array_[i + j * nx_] * exp(-r * r / R) * norm;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < nx_ * ny_; ++i)
    {
        array_[i] = larray[i];
    }
    delete larray;
}

// Output figure into file
// time : if true write time on fig
void Bitmap::output(const char *fname, bool time, double t0)
{
    char dimbuff[PHP_PATH_LENGTH];

    double RGB_array[nx_ * ny_ * 3], v, vmax = 0;
    int i, j;
    ColorRGB c_red(1, 0, 0);

    if (Max_ <= 0)
    {
        vmax = 0;
        i = j = 0;
        while (i < nx_ * ny_)
        {
            v = array_[i++];
            if (v > vmax)
            {
                vmax = v;
            }
        }
    }

    i = j = 0;
    while (i < nx_ * ny_)
    {
        v = array_[i++];
        //if(v > vmax) { vmax = v;}
        v /= (Max_ > 0) ? Max_ : vmax;
        if (v > 1.0)
        {
            v = 1.0;
        }

        RGB_array[j++] = 0; // R
        RGB_array[j++] = v; // G
        RGB_array[j++] = 0; // B
    }
    if (Max_ <= 0)
        std::cerr << "vmax=" << vmax << "\n";

    sprintf(dimbuff, "%dx%d", nx_, ny_);
    //cerr<<"New Image of size "<<dimbuff<<"\n";
    //
    //cerr<<"Read Image \n";
    Image *image;
    image = new Image(dimbuff, "white");
    image->read(nx_, ny_, "RGB", DoublePixel, RGB_array);
    //cerr<<"Write Image into "<<fname<<"\n";

    // Draw time on figure
    if (time)
    {
        string t_buff;
        stringstream ss;
        ss << (t_ - t0);
        ss >> t_buff;
        image->fontPointsize(16);
        image->font((const string) "Helvetica");
        image->strokeColor(c_red);
        image->strokeWidth(1);
        ;
        image->draw(DrawableText(6, 22, t_buff));
    }
    //image.fileName(fname);
    image->write(fname);

    delete image;
}

void Bitmap::save(const char *fname)
{
    ofstream ofs(fname, ofstream::binary);

    ofs.write((char *)&nx_, sizeof(nx_));
    ofs.write((char *)&ny_, sizeof(ny_));
    ofs.write((char *)&Lx_, sizeof(Lx_));
    ofs.write((char *)&Ly_, sizeof(Ly_));
    ofs.write((char *)&t_, sizeof(t_));
    ofs.write((char *)array_, sizeof(double) * nx_ * ny_);

    ofs.close();
}
