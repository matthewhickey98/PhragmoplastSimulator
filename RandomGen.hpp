#ifndef RandomGen_HPP
#define RandomGen_HPP
#include "php_const.hpp"
#include <cmath>

const int NTAB = 32;

class RandomGen
{
  public:
    RandomGen(double min = 0.0, double max = 1.0, long seed = -1);
    double Next();
    double next_0_1();
    double next_m1_1();
    inline double boltzman(double r)
    {
        return (-log(1 - Next0to1dbl()) / r);
    }

    double Normal(double var);
    double NormalCapped(double var, double max);
    void Normal2D(double var, std::vector<double> &v);
    void Normal2DCapped(double var, std::vector<double> &v, double max);
    void Normal3D(double var, std::vector<double> &v);
    void Normal3DCapped(double var, std::vector<double> &v, double max);

    double Trapezium(double b);
    double Next0to1dbl();
    inline double next_m1_1_dbl()
    {
        return (-1.0 + 2.0 * Next0to1dbl());
    }
    inline double next_dbl()
    {
        return (min_ + (max_ - min_) * Next0to1dbl());
    }

  private:
    long rand_seed1_;
    long rand_seed2_;
    long iy_;
    long iv_[NTAB];
    double min_;
    double max_;
    int normal_iset_;
    double normal_gset_;
    int normal_cap_iset_;
    double normal_cap_gset_;
};

class RandomGen_int : public RandomGen
{
  public:
    RandomGen_int(int imin = 0, int imax = 100, long seed = -1) : RandomGen(0.0, 1.0, seed), imin_(imin), imax_(imax){};
    int next();

  private:
    // long rand_seed1_;
    // long rand_seed2_;
    // long iy_;
    // long iv_[NTAB];
    int imin_;
    int imax_;
};

#endif
