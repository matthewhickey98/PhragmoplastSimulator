#ifndef MICROTUBLE_H
#define MICROTUBLE_H
#include "RandomGen.h"
#include "FrameName.h"
#include <iostream>
#include <fstream>
#include "php_const.h"
#include "Bitmap.h"

class Phragmoplast;

// A single microtubule
class Microtuble {
public:
    Microtuble() : bleached_x1_(-1), bleached_x2_(-1), type_(PHP_TYPE_NONE), php_(0) {};
    ~Microtuble() {};

    static void set_bounds(double xmax, double ymax, double dx0)
    {
        xmax__ = xmax;
        ymax__ = ymax;
        dx0__ = dx0;
        min_L__ = dx0__ * 0.05;
        NoMature_mode__ = false;
    }

    static void set_FRAP(double x1, double y1, double x2, double y2)
    {
        Microtuble::FRAPx1__ = x1;
        Microtuble::FRAPy1__ = y1;
        Microtuble::FRAPx2__ = x2;
        Microtuble::FRAPy2__ = y2;
    }
    static bool is_inside_FRAP(double x, double y)
    {
        return ((FRAPx1__ <= x) && (x <= FRAPx2__) &&
                (FRAPy1__ <= y) && (y <= FRAPy2__));
    }

    // return y coordinate of cross with left edge
    // Return -1 if none;
    double cross_left_FRAP()
    {
        double yl = line_y(FRAPx1__);
        //std::cerr<<" yl="<<yl<<"\n";
        if ((((y1_ <= yl + 1e-8) && (yl <= y2_ + 1e-8)) ||
                ((y2_ <= yl + 1e-8) && (yl <= y1_ + 1e-8))) &&
                (FRAPy1__ <= yl) && (yl <= FRAPy2__))
        {
            return (yl);
        }
        //std::cerr<<"NO cross_left_FRAP yl="<<yl<<"\n";
        return (-1);
    }

    // return y coordinate of cross with right edge
    // Return -1 if none;
    double cross_right_FRAP()
    {
        double yr = line_y(FRAPx2__);
        //std::cerr<<" yr="<<yr<<"\n";
        if ((((y1_ <= yr + 1e-8) && (yr <= y2_ + 1e-8)) ||
                ((y2_ <= yr + 1e-8) && (yr <= y1_ + 1e-8))) &&
                (FRAPy1__ <= yr) && (yr <= FRAPy2__))
        {
            return (yr);
        }
        //std::cerr<<"NO cross_right_FRAP yr="<<yr<<"\n";
        return (-1);
    }

    // return y coordinate of cross with bottom edge
    // Return -1 if none;
    double cross_bottom_FRAP()
    {
        double xb = line_x(FRAPy1__);
        //std::cerr<<" xb="<<xb<<"\n";
        if ((((x1_ <= xb + 1e-8) && (xb <= x2_ + 1e-8)) ||
                ((x2_ <= xb + 1e-8) && (xb <= x1_ + 1e-8))) &&
                (FRAPx1__ <= xb) && (xb <= FRAPx2__))
        {
            return (xb);
        }
        //std::cerr<<"NO cross_bottom_FRAP xb="<<xb<<"\n";
        return (-1);
    }

    // return y coordinate of cross with top edge
    // Return -1 if none;
    double cross_top_FRAP()
    {
        double xt = line_x(FRAPy2__);
        //std::cerr<<" xt="<<xt<<"\n";
        if ((((x1_ <= xt + 1e-8) && (xt <= x2_ + 1e-8)) ||
                ((x2_ <= xt + 1e-8) && (xt <= x1_ + 1e-8))) &&
                (FRAPx1__ <= xt) && (xt <= FRAPx2__))
        {
            return (xt);
        }
        //std::cerr<<"NO cross_bottom_FRAP xt="<<xt<<"\n";
        return (-1);
    }

    void init(double x1, double y1, double L, double theta, int type,
              Phragmoplast *php);
    void check_mt_pos();
    void update_bleach();
    void size_adjust(bool plus = true, bool grow = true);
    void do_event(int event);

    int select_next_event(double &dt);
    void run_until(double t_max);

    inline double get_x1()
    {
        return (x1_);
    };
    inline double get_x2()
    {
        return (x2_);
    };
    inline double get_y1()
    {
        return (y1_);
    };
    inline double get_y2()
    {
        return (y2_);
    };
    inline double L()
    {
        return ((this->state_p_ == PHP_STATE_EMPTY) ? 0 : sqrt((this->x2_ - this->x1_) * (this->x2_ - this->x1_) + (this->y2_ - this->y1_) * (this->y2_ - this->y1_)));
    };

    void draw(Bitmap &bm);
    void bleach();
    void txt_save(std::ofstream &ofs);

    inline int get_type()
    {
        return (type_);
    }
    inline void set_type(int type)
    {
        type_ = type;
    }
    inline int get_state_p()
    {
        return (state_p_);
    }
    inline int get_state_m()
    {
        return (state_m_);
    }
    inline void set_state_p(int state)
    {
        state_p_ = state;
    }
    inline void set_state_m(int state)
    {
        state_m_ = state;
    }
    inline double get_t()
    {
        return (t_);
    }
    inline void set_t(double t)
    {
        t_ = t;
    }
    double rnd(double x);

    // corrdinate of the MT
    double line_x(double y)
    {
        return (((x2_ - x1_) * y + x1_ * y2_ - x2_ * y1_) / (y2_ - y1_));
    }
    double line_y(double x)
    {
        return (((y2_ - y1_) * x + y1_ * x2_ - y2_ * x1_) / (x2_ - x1_));
    }

    static double xmax__, ymax__, dx0__, min_L__, theta_max__;
    static double FRAPx1__, FRAPy1__, FRAPx2__, FRAPy2__;
    static bool NoMature_mode__;

protected:
    double x1_;          // position of the seed
    double x2_;          // end of growing MT
    double y1_;          // position of the seed
    double y2_;          // end of growing MT
    double dx_;          // dx increment
    double dy_;          // dy increment
    double theta_;       // orientation
    double bleached_x1_; // Left coordinate bleached region
    double bleached_x2_; // Right coordinate bleached region
    double t_;
    int type_;        // TREADMILL, SEEDED (GROWING_OUTSIDE, GROWING_INSIDE)
    bool towards_DZ_; // true -> +end away from CP
    int state_p_;     // +end : EMPTY, GROWING, SHRINKING, MATURE
    int state_m_;     // -end : EMPTY, GROWING, SHRINKING, MATURE
    Phragmoplast *php_;
};

#endif
