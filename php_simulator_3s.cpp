//#include <boost/python.hpp>
#include "Phragmoplast.h"
#include <Magick++.h>

#include "RandomGen.h"
#include "FrameName.h"
#include "php_const.h"
#include "Bitmap.h"
#include "Microtuble.h"
#include "readpars.h"

const char *VERSION = "php_simulator: v1.1; 23/1/20; B. Piette";

struct Pars
{
    double Thickness;                          // Thichness of Phragmoplast (length of MT in microns)
    double Width;                              // Width of Phragmoplast (radial width) (microns)
    double dx;                                 //  Microtubule growing increment
    int N_MT;                                  // Number of microtubules
    int Npixel_X;                              // Number of pixels in horizontal direction
    int Npixel_Y;                              // Number of pixels in vertical direction
    double t_min;                              // If < 0 : relax until t=0. No data output before t=0.
    double t_max;                              // simulation duration in sec
    double FRAP_t;                             // FRAP parameters time
    double FRAP_x1, FRAP_y1, FRAP_x2, FRAP_y2; // FRAP rectangle
    double grf_dt;                             // time steps between plots, in sec
    double grf_max;                            // Max no of supperposition. Mapped to G=1
    double kymograph_amp;                      //  Kymograph Amplification
    double L0;                                 //  Initial length

    int bin_N; // Number of Bins for histogram

    double r_polym;   // growing rate in micron/s
    double r_depolym; //shrinking rate in micron/s
    double r_gs;
    double r_sg;
    double r_reseed;

    double r_pg;    // Rate pause -> growing
    double r_ps;    // Rate pause -> shrinking
    double r_gp;    // Rate growing -> pause
    double r_sp;    // Rate shrinking -> pause
    double r_ps_CP; // Rate catastrophe at Cell Plate
    double r_ps_DZ; // Rate catastrophe at Distal Zone

    double r_me_polym;   // -end Polymerisation Rate
    double r_me_depolym; // -end dePolymerisation Rate
    double r_me_gs;      // -end Catastrophe Rate
    double r_me_sg;      // -end Recovery Rate
    double r_me_pg;      // -end Rate pause -> growing
    double r_me_ps;      // -end Rate pause -> shrinking
    double r_me_gp;      // -end Rate growing -> pause
    double r_me_sp;      // -end Rate shrinking -> pause

    int NoMature_mode; // No Mature -> Growing

    double TdensCorr;
    double TaxolCorr;

    double frac_treadmill;         // fraction of treadmilling MT
    double frac_grow_out;          // fraction of out growing MT
    double frac_grow_in;           // fraction of in growing MT
    double theta_max;              // Max theta angle for MT orientation
    double frac_seed_CP;           // Fraction of CP seeded MT
    double frac_seed_distal;       // Fraction of outside edge seeded MT
    double frac_seed_middle;       // Fraction of homogeneously distributed MT.
    double frac_seed_middle_slope; // Assymetry dist. range [0-2].

    int rnd_seed; // seed for random generator (< 0 -> use time).
    char fname_prefix[PHP_PATH_LENGTH];
    char init_type[PHP_PATH_LENGTH];
    int no_grf; // Do not ouput grf files

    double dens_MT;      // Total length (in um) of MT per square um
    double macet_dz_len; // length of DZ for macet
    double macet_mz;     // macet_factor 0 -> Width-width_dz
    double macet_dz;     // macet_factor at Width
    double macet_MT_len; // Max MT lenght affected by macet

} par;

// enum ParType
// {
//   T_NONE = 0,
//   T_INT = 1,
//   T_DOUBLE = 3,
//   T_STRING = 4,
// };

// struct Parameter
// {
//   const char *name;
//   ParType type;
//   void *p;
//   const char *help;
// };

Parameter Pars[] =
{
    {"Thickness", T_DOUBLE, &par.Thickness, " "},
    {"Width", T_DOUBLE, &par.Width, " Radial width in micron"},
    {"dx", T_DOUBLE, &par.dx, " Step growth in micron"},
    {"N_MT", T_INT, &par.N_MT, " Number of microtubules"},
    {"Npixel_X", T_INT, &par.Npixel_X, " Number of pixels, horizontal"},
    {"Npixel_Y", T_INT, &par.Npixel_Y, " Number of pixels, vertical"},
    {"t_min", T_DOUBLE, &par.t_min, " "},
    {"t_max", T_DOUBLE, &par.t_max, " "},
    {"FRAP_t", T_DOUBLE, &par.FRAP_t, " "},
    {"FRAP_x1", T_DOUBLE, &par.FRAP_x1, " "},
    {"FRAP_y1", T_DOUBLE, &par.FRAP_y1, " "},
    {"FRAP_x2", T_DOUBLE, &par.FRAP_x2, " "},
    {"FRAP_y2", T_DOUBLE, &par.FRAP_y2, " "},
    {"grf_dt", T_DOUBLE, &par.grf_dt, " "},
    {"grf_max", T_DOUBLE, &par.grf_max, "Max no of supperposition. Mapped to G=1"},
    {"kymograph_amp", T_DOUBLE, &par.kymograph_amp, "kymograph amplitude."},
    {"L0", T_DOUBLE, &par.L0, " "},
    {"bin_N", T_INT, &par.bin_N, " Number of histogram bins"},
    {"r_polym", T_DOUBLE, &par.r_polym, "growing rate in micron/s"},
    {"r_depolym", T_DOUBLE, &par.r_depolym, "shrinking rate in micron/s"},
    {"r_gs", T_DOUBLE, &par.r_gs, " catastrophe rate "},
    {"r_sg", T_DOUBLE, &par.r_sg, " Recovery rate"},
    {"r_reseed", T_DOUBLE, &par.r_reseed, " Reseed rate"},
    {"r_pg", T_DOUBLE, &par.r_pg, " Rate pause -> growing"},
    {"r_ps", T_DOUBLE, &par.r_ps, " Rate pause -> shrinking"},
    {"r_gp", T_DOUBLE, &par.r_gp, " Rate growing -> pause"},
    {"r_sp", T_DOUBLE, &par.r_sp, " Rate shrinking -> pause"},
    {"r_ps_CP", T_DOUBLE, &par.r_ps_CP, " Rate catastrophe at Cell Plate"},
    {"r_ps_DZ", T_DOUBLE, &par.r_ps_DZ, " Rate catastrophe at Distal Zone"},

    {"r_me_polym", T_DOUBLE, &par.r_me_polym, " -end polym rate"},
    {"r_me_depolym", T_DOUBLE, &par.r_me_depolym, " -end depolym rate"},
    {"r_me_gs", T_DOUBLE, &par.r_me_gs, " -end catastrophe rate"},
    {"r_me_sg", T_DOUBLE, &par.r_me_sg, " -end rescue rate"},
    {"r_me_pg", T_DOUBLE, &par.r_me_pg, " -end Rate pause -> growing"},
    {"r_me_ps", T_DOUBLE, &par.r_me_ps, " -end Rate pause -> shrinking"},
    {"r_me_gp", T_DOUBLE, &par.r_me_gp, " -end Rate growing -> pause"},
    {"r_me_sp", T_DOUBLE, &par.r_me_sp, " -end Rate shrinking -> pause"},

    // ASADA EXPERIMENT
    {"NoMature_mode", T_INT, &par.NoMature_mode, " No mature state"},
    {"TdensCorr", T_DOUBLE, &par.TdensCorr, " "},
    {"TaxolCorr", T_DOUBLE, &par.TaxolCorr, " "},

    {"frac_treadmill", T_DOUBLE, &par.frac_treadmill, "Fraction MT treadmilling "},
    {"frac_grow_out", T_DOUBLE, &par.frac_grow_out, "Fraction MT groing in "},
    {"frac_grow_in", T_DOUBLE, &par.frac_grow_in, "Fraction MT growing out "},
    {"theta_max", T_DOUBLE, &par.theta_max, " Max growing angle for MT"},
    {"frac_seed_CP", T_DOUBLE, &par.frac_seed_CP, " Fraction seeds on Cell Plate "},
    {
        "frac_seed_distal", T_DOUBLE, &par.frac_seed_distal,
        "Fraction seeds on distal zone "
    },
    {
        "frac_seed_middle", T_DOUBLE, &par.frac_seed_middle,
        "Fraction seeds in mid zone "
    },
    {
        "frac_seed_middle_slope", T_DOUBLE, &par.frac_seed_middle_slope,
        "Fraction seeds in mid zone slope "
    },

    {
        "rnd_seed", T_INT, &par.rnd_seed,
        "seed for random generator (< 0 -> use time)."
    },
    {"fname_prefix", T_STRING, par.fname_prefix, " graphic file prefix"},
    {"init_type", T_STRING, par.init_type, " RND, FIXED_L, DEBUG"},
    {"no_grf", T_INT, &par.no_grf, " do not output grf files"},

    {"dens_MT", T_DOUBLE, &par.dens_MT, " Length density (um) of MT per um^2"},
    {"macet_dz_len", T_DOUBLE, &par.macet_dz_len, " Length of DZ for macet"},
    {"macet_mz", T_DOUBLE, &par.macet_mz, " macet_factor 0 -> Width-width_dz"},
    {"macet_dz", T_DOUBLE, &par.macet_dz, " macet_factor at Width"},
    {"macet_MT_len", T_DOUBLE, &par.macet_MT_len, " Max MT lenght affected by macet"},

    {0, T_NONE, 0},
};

const char *help[] =
{
    "php_simulator PAR_FILE",
    0
};

const char *init_types[] =
{
    "RND",
    "FIXED_L",
    "DEBUG",
    0,
};
const int INIT_TYPE_RND = 0;
const int INIT_TYPE_FIXED_L = 1;
const int INIT_TYPE_DEBUG = 2;


// BOOST_PYTHON_MODULE(simulator)
// {
//     def("runSim", simMain);
// }



int main(int argc, char **argv)
{
    int init_type;
    Magick::InitializeMagick(NULL);
    try
    {
        read_pars(argv[1], Pars);
    }
    catch(const std::exception& e)
    {
        cout << "Exception Caught: " << e.what() << endl;
        exit(1);
    }

    if ((init_type = string_index(init_types, par.init_type)) < 0)
    {
        std::cerr << "Invalid init_type :" << par.init_type << "\n";
        exit(0);
    }

    Phragmoplast php(par.N_MT, par.Thickness, par.Width, par.dx,
                     par.Npixel_X, par.Npixel_Y, par.grf_dt,
                     par.r_polym * par.TdensCorr, par.r_depolym,
                     par.r_gs / (par.TdensCorr * par.TaxolCorr),
                     par.r_sg * par.TdensCorr * par.TaxolCorr,
                     par.r_reseed,
                     par.r_pg * par.TdensCorr * par.TaxolCorr,
                     par.r_ps / (par.TdensCorr * par.TaxolCorr),
                     par.r_gp / (par.TdensCorr * par.TaxolCorr),
                     par.r_sp * par.TdensCorr * par.TaxolCorr,
                     par.r_ps_CP / (par.TdensCorr * par.TaxolCorr),
                     par.r_ps_DZ / (par.TdensCorr * par.TaxolCorr),
                     par.r_me_polym * par.TdensCorr * par.TaxolCorr,
                     par.r_me_depolym,
                     par.r_me_gs / (par.TdensCorr * par.TaxolCorr),
                     par.r_me_sg * par.TdensCorr * par.TaxolCorr,
                     par.r_me_pg * par.TdensCorr * par.TaxolCorr,
                     par.r_me_ps / (par.TdensCorr * par.TaxolCorr),
                     par.r_me_gp / (par.TdensCorr * par.TaxolCorr),
                     par.r_me_sp * par.TdensCorr * par.TaxolCorr,
                     par.fname_prefix, par.grf_max, par.kymograph_amp,
                     par.dens_MT, par.macet_dz_len, par.macet_mz, par.macet_dz,
                     par.macet_MT_len);

    php.set_bin_N(par.bin_N);

    if (par.NoMature_mode)
    {
        Microtuble::NoMature_mode__ = true;
    }

    switch (init_type)
    {
        case INIT_TYPE_RND:
            php.init(par.frac_treadmill, par.frac_grow_out, par.frac_grow_in,
                     par.theta_max, par.frac_seed_CP, par.frac_seed_distal,
                     par.frac_seed_middle, par.frac_seed_middle_slope,
                     par.L0); // NEW added extra par here
            break;

        case INIT_TYPE_FIXED_L:
            php.init_fixed_L(par.L0, par.theta_max);
            break;
        default:
            exit(0);
    }
    php.set_FRAP(par.FRAP_t, par.FRAP_x1, par.FRAP_y1, par.FRAP_x2, par.FRAP_y2);
    if (par.no_grf)
    {
        php.set_no_grf();
    }

    php.run_sync(par.t_min, par.t_max, par.rnd_seed);

}
