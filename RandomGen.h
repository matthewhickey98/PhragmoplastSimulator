#ifndef RandomGen_H
#define RandomGen_H
#include <cmath>

const int NTAB = 32;

class RandomGen {
public:
    RandomGen(double min = 0.0, double max = 1.0, long seed = -1);
    double next();
    double next_0_1();
    double next_m1_1();
    inline double boltzman(double r)
    {
        return (-log(1 - next_0_1_dbl()) / r);
    }

    double normal(double var);
    double normal_capped(double var, double max);
    void normal_2d(double var, double *v);
    void normal_2d_capped(double var, double *v, double max);
    void normal_3d(double var, double *v);
    void normal_3d_capped(double var, double *v, double max);

    double trapezium(double b);

    double next_0_1_dbl();
    inline double next_m1_1_dbl()
    {
        return (-1.0 + 2.0 * next_0_1_dbl());
    }
    inline double next_dbl()
    {
        return (min_ + (max_ - min_) * next_0_1_dbl());
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

class RandomGen_int : public RandomGen {
public:
    RandomGen_int(int imin = 0, int imax = 100, long seed = -1) : RandomGen(0.0, 1.0, seed), imin_(imin), imax_(imax) {};
    int next();

private:
    long rand_seed1_;
    long rand_seed2_;
    long iy_;
    long iv_[NTAB];
    int imin_;
    int imax_;
};

#endif
