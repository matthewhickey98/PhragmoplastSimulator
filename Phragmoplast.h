#ifndef PHP_SIM_H
#define PHP_SIM_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <Magick++.h>
#include <vector>

#include "RandomGen.h"
#include "FrameName.h"
#include "php_const.h"
#include "Bitmap.h"
#include "Microtuble.h"
#include "readpars.h"

// Phragmoplast made out of MT
class Phragmoplast {
public:
    Phragmoplast(int N_MT, double Thickness, double Width, double dx,
                 int Npixel_X, int Npixel_Y, double grf_dt,
                 double r_polym, double r_depolym,
                 double r_gs,
                 double r_sg, double r_reseed,
                 double r_pg, double r_ps, double r_gp, double r_sp,
                 double r_ps_CP, double r_ps_DZ,
                 double r_me_polym, double r_me_depolym, double r_me_gs,
                 double r_me_sg,
                 double r_me_pg, double r_me_ps, double r_me_gp, double r_me_sp,
                 const char *figname, double grf_max,
                 double kymograph_amp,
                 double dens_MT,
                 double macet_dz_len, double macet_mz, double macet_dz,
                 double macet_MT_len);
    ~Phragmoplast()
    {
        delete[] mt_array_;
        delete kymograph_bm_;
    }

    void init(double frac_treadmill, double frac_grow_out,
              double frac_grow_in, double theta_max,
              double frac_seed_CP, double frac_seed_distal,
              double frac_seed_middle, double frac_seed_middle_slope,
              double L0); // NEW : added L0 as parameter here);
    void init_fixed_L(double L0, double theta_max);

    void set_FRAP(double t, double x1, double y1, double x2, double y2)
    {
        FRAP_t_ = t;
        FRAP_x1_ = x1;
        FRAP_y1_ = y1;
        FRAP_x2_ = x2;
        FRAP_y2_ = y2;
        Microtuble::set_FRAP(x1, y1, x2, y2);
    }

    void run_sync(double t_min, double t_max, long seed = -1);
    void evolve_sync_all_MT_until(double tmax);

    void bleach();
    double total_MT_length(double &average_L);

    double rnd(double xmax)
    {
        return (rnd_.next_0_1_dbl() * xmax);
    }

    double random_t(double rate)
    {
        double r = rnd_.next_0_1();
        //std::cerr<<"rate="<<rate<<"  r="<<r<<"  t="<<-log(r)/rate<<"\n";
        //  if(r < 0 || r >= 1) { std::cerr <<"ERROR: rnd = "<<r<<"\n"; exit(0); }
        return (-log(r) / rate);
    };

    void make_fig(double t);
    void make_histogram(double t);
    void save_list(double t);
    void add_kymograph_line(Bitmap &bm);
    void make_kymograph();

    inline double get_Thickness()
    {
        return (Thickness_);
    }
    inline double get_Width()
    {
        return (Width_);
    }
    inline double get_theta_max()
    {
        return (theta_max_);
    }

    inline double get_r_polym()
    {
        return (r_polym_);
    }
    inline double get_r_depolym()
    {
        return (r_depolym_);
    }
    inline double get_r_gs()
    {
        return (r_gs_);
    }
    inline double get_r_sg()
    {
        return (r_sg_);
    }
    inline double get_r_reseed()
    {
        return (r_reseed_);
    }
    inline double get_r_pg()
    {
        return (r_pg_);
    }
    inline double get_r_ps()
    {
        return (r_ps_);
    }
    inline double get_r_gp()
    {
        return (r_gp_);
    }
    inline double get_r_sp()
    {
        return (r_sp_);
    }
    inline double get_r_ps_CP()
    {
        return (r_ps_CP_);
    }
    inline double get_r_ps_DZ()
    {
        return (r_ps_DZ_);
    }

    inline double get_r_me_polym()
    {
        return (r_me_polym_);
    }
    inline double get_r_me_depolym()
    {
        return (r_me_depolym_);
    }
    inline double get_r_me_gs()
    {
        return (r_me_gs_);
    }
    inline double get_r_me_sg()
    {
        return (r_me_sg_);
    }
    inline double get_r_me_pg()
    {
        return (r_me_pg_);
    }
    inline double get_r_me_ps()
    {
        return (r_me_ps_);
    }
    inline double get_r_me_gp()
    {
        return (r_me_gp_);
    }
    inline double get_r_me_sp()
    {
        return (r_me_sp_);
    }

    inline void set_no_grf()
    {
        no_grf_ = true;
    }
    inline void set_bin_N(int bin_N)
    {
        bin_N_ = bin_N;
    }

    // MT dependant rates
    inline double R_polym()
    {
        return (r_polym_ * rel_MT_);
    }
    double R_gs(double L_MT, double x);
    //inline double R_gs_mature() { return(r_gs_mature_/(rel_MT_+1e-6)); }
    //inline double R_gs_MZ() { return(r_gs_MZ_/(rel_MT_+1e-6)); }
    inline double R_r_reseed()
    {
        return (r_reseed_ * (rel_MT_));
    }

    inline double R_sg()
    {
        return (r_sg_ * rel_MT_);
    }
    inline double R_pg()
    {
        return (r_pg_ * rel_MT_);
    }
    inline double R_ps()
    {
        return (r_ps_ / (rel_MT_ + 1e-6));
    }
    double R_gp(double L_MT, double x);
    inline double R_sp()
    {
        return (r_sp_ * rel_MT_);
    }
    inline double R_ps_CP()
    {
        return (r_ps_CP_ / (rel_MT_ + 1e-6));
    }
    inline double R_ps_DZ()
    {
        return (r_ps_DZ_ / (rel_MT_ + 1e-6 ));
    }

    inline double R_me_polym()
    {
        return (r_me_polym_ * rel_MT_);
    }
    double R_me_gs(double L_MT, double x);
    inline double R_me_sg()
    {
        return (r_me_sg_ * rel_MT_);
    }
    inline double R_me_pg()
    {
        return (r_me_pg_ * rel_MT_);
    }
    inline double R_me_ps()
    {
        return (r_me_ps_ / (rel_MT_));
    }
    double R_me_gp(double L_MT, double x);
    inline double R_me_sp()
    {
        return (r_me_sp_ * rel_MT_);
    }

protected:
    int N_MT_;         // Number of microtubules
    double Thickness_; // Thickness of phragmoplast (length of MT in microns)
    double Width_;     // Width of phragmoplast (radial width) (microns)
    double dx_;        // MT step growth (micron, use 100 nm?)
    double theta_max_;

    double t_;

    int Npixel_X_;  // Number pixels horizontaly
    int Npixel_Y_;  // Number pixels verticaly
    double grf_dt_; // time step between figures
    double MaxSupperposition_;

    Microtuble *mt_array_;

    double r_polym_;     // micron/s
    double r_depolym_;   // shrinking rate in micron/s
    double r_gs_;        // catastrophe rate in 1/s
    double r_gs_mature_; // catastrophe rate when fully grown in 1/s
    double gs_mature_len_;
    double DZ_edge_corr_;
    double r_gs_MZ_;
    double gs_MZ_len_;
    double MZ_edge_corr_;
    double r_sg_;     // recovery rate in 1/s
    double r_reseed_; // reseeding rate in 1/s

    double r_pg_;    // Rate pause -> growing
    double r_ps_;    // Rate pause -> shrinking
    double r_gp_;    // Rate growing -> pause
    double r_sp_;    // Rate shrinking -> pause
    double r_ps_CP_; // Rate shrinking -> pause at CP
    double r_ps_DZ_; // Rate shrinking -> pause at DZ

    // Rates for -end
    double r_me_polym_;   // Polymerisation rate
    double r_me_depolym_; // Depolymerisation rate
    double r_me_gs_;      // Catastrophe rate
    double r_me_sg_;      // Rescure rate
    double r_me_pg_;      // Rate pause -> growing
    double r_me_ps_;      // Rate pause -> shrinking
    double r_me_gp_;      // Rate growing -> pause
    double r_me_sp_;      // Rate shrinking -> pause

    double FRAP_t_; // Frap parameters
    double FRAP_x1_, FRAP_y1_, FRAP_x2_, FRAP_y2_;

    char fig_name_[PHP_PATH_LENGTH + 30];
    char kymograph_name_[PHP_PATH_LENGTH + 60];
    char luminosity_name_[PHP_PATH_LENGTH + 60];
    char histogram_name_[PHP_PATH_LENGTH + 60];
    FrameName fname_;
    Bitmap *kymograph_bm_;
    int index_kymograph_;
    int kymograph_lines_;
    double kymograph_amp_;
    bool no_grf_;
    int bin_N_;

public:
    RandomGen rnd_;
    double dens_MT_; // length (um) density of MT per um^2
    double macet_dz_len_;
    double macet_mz_;
    double macet_dz_;
    double macet_MT_len_;
    double tot_MT_0_; // Initial total length of MT
    double tot_MT_;   // Total length of MT
    double rel_MT_;   // tot_MT_/tot_MT_0_
    
};
#endif
