#include "RandomGen.h"
#include <ctime>
#include <cmath>

const int IM1 = 2147483563;
const int IM2 = 2147483399;
const int IA1 = 40014;
const int IA2 = 40692;
const int IQ1 = 53668;
const int IQ2 = 53774;
const int IR1 = 12211;
const int IR2 = 3791;
const int IMM1 = (IM1 - 1);
const double EPS = 1.2e-7;
//const double EPS = 1e-16;
const double AM = (1.0 / IM1);
const double NDIV = (1 + (double)IMM1 / (double)NTAB);
const double RNMAX = (1.0 - EPS);

/******************************************/
/* Initialise the random number generator */
/* min : lower bound for random number    */
/* max : upper bound for random number    */
/* seed : initial seed.                   */
/*        -1 -> use current time          */
/******************************************/
RandomGen::RandomGen(double min, double max, long seed) : rand_seed1_(0), rand_seed2_(123456789), min_(min), max_(max), normal_iset_(0), normal_cap_iset_(0)
{
    int j;
    long k;

    /* by default we take the current time as the seed */
    if (seed < 1)
        rand_seed1_ = (long)time(0);
    else
        rand_seed1_ = seed;

    rand_seed1_ = rand_seed1_ & 0xffffffff;
    rand_seed2_ = rand_seed1_;

    for (j = NTAB + 7; j >= 0; j--)
    {
        k = rand_seed1_ / IQ1;
        rand_seed1_ = IA1 * (rand_seed1_ - k * IQ1) - k * IR1;
        if (rand_seed1_ < 0)
            rand_seed1_ += IM1;
        if (j < NTAB)
            iv_[j] = rand_seed1_;
    }
    iy_ = iv_[0];
}

/*************************************************/
/* Return the next number in the range ]min,max[ */
/*************************************************/
double RandomGen::next()
{
    int j;
    long k;
    double temp;

    k = rand_seed1_ / IQ1;
    rand_seed1_ = IA1 * (rand_seed1_ - k * IQ1) - k * IR1;
    if (rand_seed1_ < 0)
        rand_seed1_ += IM1;
    k = rand_seed2_ / IQ2;
    rand_seed2_ = IA2 * (rand_seed2_ - k * IQ2) - k * IR2;
    if (rand_seed2_ < 0)
        rand_seed2_ += IM2;

    j = ((int)(iy_ / NDIV)) % NTAB;
    iy_ = iv_[j] - rand_seed2_;
    iv_[j] = rand_seed1_;
    if (iy_ < 1)
        iy_ += IMM1;
    if ((temp = AM * iy_) > RNMAX)
        temp = RNMAX;
    if (temp <= 0)
        temp = EPS;
    // temp is in the range ]0,1[

    return (min_ + (max_ - min_) * temp);
}

/*********************************************/
/* Return the next number in the range ]0,1[ */
/*********************************************/
double RandomGen::next_0_1()
{
    int j;
    long k;
    double temp;

    k = rand_seed1_ / IQ1;
    rand_seed1_ = IA1 * (rand_seed1_ - k * IQ1) - k * IR1;
    if (rand_seed1_ < 0)
        rand_seed1_ += IM1;
    k = rand_seed2_ / IQ2;
    rand_seed2_ = IA2 * (rand_seed2_ - k * IQ2) - k * IR2;
    if (rand_seed2_ < 0)
        rand_seed2_ += IM2;

    j = ((int)(iy_ / NDIV)) % NTAB;
    iy_ = iv_[j] - rand_seed2_;
    iv_[j] = rand_seed1_;
    if (iy_ < 1)
        iy_ += IMM1;
    if ((temp = AM * iy_) > RNMAX)
        temp = RNMAX;
    if (temp <= 0)
        temp = EPS;
    // temp is in the range ]0,1[

    return (temp);
}

/*************************************************/
/* Return the next number in the range ]-1,1[    */
/*************************************************/
double RandomGen::next_m1_1()
{
    int j;
    long k;
    double temp;

    k = rand_seed1_ / IQ1;
    rand_seed1_ = IA1 * (rand_seed1_ - k * IQ1) - k * IR1;
    if (rand_seed1_ < 0)
        rand_seed1_ += IM1;
    k = rand_seed2_ / IQ2;
    rand_seed2_ = IA2 * (rand_seed2_ - k * IQ2) - k * IR2;
    if (rand_seed2_ < 0)
        rand_seed2_ += IM2;

    j = ((int)(iy_ / NDIV)) % NTAB;
    iy_ = iv_[j] - rand_seed2_;
    iv_[j] = rand_seed1_;
    if (iy_ < 1)
        iy_ += IMM1;
    if ((temp = AM * iy_) > RNMAX)
        temp = RNMAX;
    if (temp <= 0)
        temp = EPS;
    // temp is in the range ]0,1[

    return (-1 + 2 * temp);
}

/******************************************/
/* Return double with increament 1e-14 dx */
/* in range  ]1e-14,1-1e-14[ */
/******************************************/
const double EPS_DBL = 1e-14;
const double RNMAX_DBL = 1 - EPS_DBL;
double RandomGen::next_0_1_dbl()
{
    double d2 = next_m1_1();
    d2 = (next_0_1() - EPS) * (1 - 2 * EPS_DBL) / (1 - 2 * EPS) + EPS_DBL + d2 * 1.2e-7;
    if (d2 < EPS_DBL)
    {
        d2 += 1 - 2 * EPS_DBL;
    }
    else if (d2 > RNMAX_DBL)
    {
        d2 -= 1 - 2 * EPS_DBL;
    }

    return (d2);
}

/*****************************************************/
/* Normal distribution of zero mean and variance var */
/*****************************************************/
double RandomGen::normal(double var)
{
    double fac, r, v1, v2;

    if (normal_iset_ == 0)
    {
        do
        {
            v1 = next_m1_1();
            v2 = next_m1_1();
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = var * sqrt(-2.0 * log(r) / r);

        normal_gset_ = v1 * fac;
        normal_iset_ = 1;
        return (v2 * fac);
    }
    else
    {
        normal_iset_ = 0;
        return (normal_gset_);
    }
}

/*****************************************************/
/* Normal distribution of zero mean and variance var */
/*****************************************************/
double RandomGen::normal_capped(double var, double max)
{
    double fac, r, v1, v2, v;

    if (normal_cap_iset_ == 0)
    {
        do
        {
            v1 = next_m1_1();
            v2 = next_m1_1();
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = var * sqrt(-2.0 * log(r) / r);

        normal_cap_gset_ = v1 * fac;
        normal_cap_iset_ = 1;
        v = v2 * fac;
    }
    else
    {
        normal_cap_iset_ = 0;
        v = normal_cap_gset_;
    }
    if (v > max)
        v = max;
    return (v);
}

/******************************/
/* Return normal 2 dim vector */
/******************************/
void RandomGen::normal_2d(double var, double *v)
{
    double r, phi;

    r = normal(var);
    phi = next_0_1() * 2 * M_PI;
    v[0] = r * cos(phi);
    v[1] = r * sin(phi);
}

/******************************/
/* Return normal 2 dim vector */
/******************************/
void RandomGen::normal_2d_capped(double var, double *v, double max)
{
    double r, phi;

    r = normal_capped(var, max);
    phi = next_0_1() * 2 * M_PI;
    v[0] = r * cos(phi);
    v[1] = r * sin(phi);
}

/******************************/
/* Return normal 3 dim vector */
/******************************/
void RandomGen::normal_3d(double var, double *v)
{
    double r, phi, theta;

    r = normal(var);
    phi = next_0_1() * 2 * M_PI;
    theta = next_0_1() * M_PI;
    v[0] = r * sin(theta) * cos(phi);
    v[1] = r * sin(theta) * sin(phi);
    v[2] = r * cos(theta);
}

/******************************/
/* Return normal 3 dim vector */
/******************************/
void RandomGen::normal_3d_capped(double var, double *v, double max)
{
    double r, phi, theta;

    r = normal_capped(var, max);
    phi = next_0_1() * 2 * M_PI;
    theta = next_0_1() * M_PI;
    v[0] = r * sin(theta) * cos(phi);
    v[1] = r * sin(theta) * sin(phi);
    v[2] = r * cos(theta);
}

/********************************************/
/* Trapezium distribution                   */
/*                                          */
/*  |                ----|                  */
/*  |            ----    |                  */
/*  |        ----        |                  */
/*  |    ----            |                  */
/*  |----                |                  */
/*  |                    |                  */
/*  |                    |                  */
/* Slope : b   0 < b <= 2                   */
/* p(y)=1-b/2+b*y;                          */
/*   => x = int p(y) dy = (1-b/2)*y+b*y^2/2 */
/*   => y = (b/2-1+sqrt((1-b/2)^2+2bx)))/b  */
/********************************************/
double RandomGen::trapezium(double b)
{
    if ((-1e-7 < b) && (b < 1e-7))
    {
        return (next());
    } // slope too smal
    return ((b * 0.5 - 1 + sqrt((b * 0.5 - 1) * (b * 0.5 - 1) + 2 * b * next())) / b);
}

/***************************************/
/* Return a random number in the range */
/* [imin, imax ]                       */
/***************************************/
int RandomGen_int::next()
{
    return (imin_ + (int)floor(RandomGen::next() * (imax_ - imin_ + 1)));
}
