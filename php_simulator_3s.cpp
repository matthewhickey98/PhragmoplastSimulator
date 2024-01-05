//#include <boost/python.hpp>
#include "Bitmap.hpp"
#include "FrameName.hpp"
#include "Magick++.h"
#include "Microtubule.hpp"
#include "Phragmoplast.hpp"
#include "RandomGen.hpp"
#include "php_const.hpp"
#include "readpars.hpp"

const char *VERSION = "php_simulator: v2; 10/25/2021; M.Hickey, B. Piette";

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
    double r_depolym; // shrinking rate in micron/s
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
    int SimIndex;
    double SeedParameterA;
    double SeedParameterB;
    int distType; // Type 0 -> Cauchy | Type 1 -> Normal
    double CP_X_MAX;
    double DZ_X_MIN;
    double DZ_SEED_MIN;
    double reseedSlope;
    double reseedB;
    int outputJPEG;
    double blurRadiusJPEG;
    int outputBMA;
    int writeTime;

} par;

sParameter Pars[] = {
    {"fname_prefix", T_STRING, par.fname_prefix, "graphic file prefix"},
    {"DistType", T_INT, &par.distType, "Distribution type. 0 -> Cauchy | 1 -> Normal"},
    {"NoMature_mode", T_INT, &par.NoMature_mode, " No mature state"},
    {"Thickness", T_DOUBLE, &par.Thickness, "Thickness in microns"},
    {"Width", T_DOUBLE, &par.Width, " Radial width in microns"},
    {"dx", T_DOUBLE, &par.dx, " Step growth in micron"},
    {"N_MT", T_INT, &par.N_MT, " Number of microtubules"},
    {"Npixel_X", T_INT, &par.Npixel_X, " Number of pixels, horizontal"},
    {"Npixel_Y", T_INT, &par.Npixel_Y, " Number of pixels, vertical"},
    {"t_min", T_DOUBLE, &par.t_min, "Minimum simulation time"},
    {"t_max", T_DOUBLE, &par.t_max, "Maximum simulation time"},
    {"FRAP_t", T_DOUBLE, &par.FRAP_t, "FRAP Time"},
    {"FRAP_x1", T_DOUBLE, &par.FRAP_x1, "FRAP X1 Coordinate"},
    {"FRAP_y1", T_DOUBLE, &par.FRAP_y1, "FRAP Y1 Coordinate"},
    {"FRAP_x2", T_DOUBLE, &par.FRAP_x2, "FRAP X2 Coordinate"},
    {"FRAP_y2", T_DOUBLE, &par.FRAP_y2, "FRAP Y2 Coordinate"},
    {"bin_N", T_INT, &par.bin_N, " Number of histogram bins"},
    {"outputJPEG", T_INT, &par.outputJPEG, "1: to generate JPEG images"},
    {"outputBMA", T_INT, &par.no_grf, "1: to generate BMA images"},
    {"blurRadiusJPEG", T_DOUBLE, &par.blurRadiusJPEG, "Blur radius for JPEG images"},
    {"init_type", T_STRING, par.init_type, "RND, FIXED_L, DEBUG"},
    {"writeImageTime", T_INT, &par.writeTime, "Write time on figures"},
    {"rnd_seed", T_INT, &par.rnd_seed, "seed for random generator (< 0 -> use time)."},
    {"SimIndex", T_INT, &par.SimIndex, "Gloabl Index for MATLAB Sims"},

    {"r_polym", T_DOUBLE, &par.r_polym, "growing rate in micron/s"},
    {"r_depolym", T_DOUBLE, &par.r_depolym, "shrinking rate in micron/s"},
    {"r_sg", T_DOUBLE, &par.r_sg, "Recovery rate"},
    {"r_gs", T_DOUBLE, &par.r_gs, "cCtastrophe rate "},
    {"r_pg", T_DOUBLE, &par.r_pg, "Rate pause -> growing"},
    {"r_ps", T_DOUBLE, &par.r_ps, "Rate pause -> shrinking"},
    {"r_gp", T_DOUBLE, &par.r_gp, "Rate growing -> pause"},
    {"r_sp", T_DOUBLE, &par.r_sp, "Rate shrinking -> pause"},
    {"r_ps_CP", T_DOUBLE, &par.r_ps_CP, "Rate catastrophe at Cell Plate"},
    {"r_ps_DZ", T_DOUBLE, &par.r_ps_DZ, "Rate catastrophe at Distal Zone"},

    {"dens_MT", T_DOUBLE, &par.dens_MT, "Length density (um) of MT per um^2"},
    {"L0", T_DOUBLE, &par.L0, "Initial Length/Maximum intial length"},
    {"r_reseed", T_DOUBLE, &par.r_reseed, " Reseed rate"},
    {"frac_seed_middle_slope", T_DOUBLE, &par.frac_seed_middle_slope, "Fraction seeds in mid zone slope "},

    {"frac_grow_out", T_DOUBLE, &par.frac_grow_out, "Fraction MT groing in "},
    {"frac_grow_in", T_DOUBLE, &par.frac_grow_in, "Fraction MT growing out "},
    {"frac_seed_distal", T_DOUBLE, &par.frac_seed_distal, "Fraction seeds on distal zone "},
    {"frac_seed_middle", T_DOUBLE, &par.frac_seed_middle, "Fraction seeds in mid zone "},
    {"frac_seed_CP", T_DOUBLE, &par.frac_seed_CP, " Fraction seeds on Cell Plate "},
    {"frac_treadmill", T_DOUBLE, &par.frac_treadmill, "Fraction MT treadmilling "},
    {"theta_max", T_DOUBLE, &par.theta_max, " Max growing angle for MT"},
    {"CP_X_MAX", T_DOUBLE, &par.CP_X_MAX, "Cell plate x boundary (end of region)"},
    {"DZ_X_MIN", T_DOUBLE, &par.DZ_X_MIN, "Distal x min (start of region) for r_ps_DZ"},
    {"DZ_SEED_MIN", T_DOUBLE, &par.DZ_SEED_MIN, "Min x for seeding in DZ zone"},
    {"grf_dt", T_DOUBLE, &par.grf_dt, "Time step between images"},
    {"grf_max", T_DOUBLE, &par.grf_max, "Max no of supperposition. Mapped to G=1"},
    {"kymograph_amp", T_DOUBLE, &par.kymograph_amp, "kymograph amplitude."},
    {"r_me_polym", T_DOUBLE, &par.r_me_polym, "-end polym rate"},
    {"r_me_depolym", T_DOUBLE, &par.r_me_depolym, "-end depolym rate"},
    {"r_me_gs", T_DOUBLE, &par.r_me_gs, "-end catastrophe rate"},
    {"r_me_sg", T_DOUBLE, &par.r_me_sg, "-end rescue rate"},
    {"r_me_pg", T_DOUBLE, &par.r_me_pg, "-end Rate pause -> growing"},
    {"r_me_ps", T_DOUBLE, &par.r_me_ps, "-end Rate pause -> shrinking"},
    {"r_me_gp", T_DOUBLE, &par.r_me_gp, "-end Rate growing -> pause"},
    {"r_me_sp", T_DOUBLE, &par.r_me_sp, "-end Rate shrinking -> pause"},

    // ASADA EXPERIMENT
    {"TaxolCorr", T_DOUBLE, &par.TaxolCorr, "Taxol correction factor"},
    {"TdensCorr", T_DOUBLE, &par.TdensCorr, "Density correction factor"},
    {"SeedParameterA", T_DOUBLE, &par.SeedParameterA, "Parameter A for seeding distribution"},
    {"SeedParameterB", T_DOUBLE, &par.SeedParameterB, "Parameter A for seeding distribution"},

    {"reseed_slope", T_DOUBLE, &par.reseedSlope, "Slope of reseed gradient"},
    {"reseed_b", T_DOUBLE, &par.reseedB, "Y intercept of reseed gradient"},

    {0, T_NONE, 0},
};

const char *help[] = {"php_simulator PAR_FILE", 0};

const char *init_types[] = {
    "RND",
    "FIXED_L",
    "DEBUG",
    0,
};

const int INIT_TYPE_RND = 0;
const int INIT_TYPE_FIXED_L = 1;
const int INIT_TYPE_DEBUG = 2;

int main(int argc, char **argv)
{
    int init_type;

    // MatFile mFile;
    // mFile.Test();

    Magick::InitializeMagick(NULL);
    try
    {
        ReadPars(argv[1], Pars);
    }
    catch (const std::exception &e)
    {
        cout << "Exception Caught: " << e.what() << endl;
        exit(1);
    }

    if ((init_type = StringIndex(init_types, par.init_type)) < 0)
    {
        std::cout << "Invalid init_type :" << par.init_type << "\n";
        exit(0);
    }

    Phragmoplast php(
        par.N_MT, par.Thickness, par.Width, par.dx, par.Npixel_X, par.Npixel_Y, par.grf_dt, par.r_polym * par.TdensCorr,
        par.r_depolym, par.r_gs / (par.TdensCorr * par.TaxolCorr), par.r_sg * par.TdensCorr * par.TaxolCorr,
        par.r_reseed, par.r_pg * par.TdensCorr * par.TaxolCorr, par.r_ps / (par.TdensCorr * par.TaxolCorr),
        par.r_gp / (par.TdensCorr * par.TaxolCorr), par.r_sp * par.TdensCorr * par.TaxolCorr,
        par.r_ps_CP / (par.TdensCorr * par.TaxolCorr), par.r_ps_DZ / (par.TdensCorr * par.TaxolCorr),
        par.r_me_polym * par.TdensCorr * par.TaxolCorr, par.r_me_depolym, par.r_me_gs / (par.TdensCorr * par.TaxolCorr),
        par.r_me_sg * par.TdensCorr * par.TaxolCorr, par.r_me_pg * par.TdensCorr * par.TaxolCorr,
        par.r_me_ps / (par.TdensCorr * par.TaxolCorr), par.r_me_gp / (par.TdensCorr * par.TaxolCorr),
        par.r_me_sp * par.TdensCorr * par.TaxolCorr, par.fname_prefix, par.grf_max, par.kymograph_amp, par.dens_MT,
        par.SimIndex, par.t_max, par.SeedParameterA,
        par.SeedParameterB, par.distType, par.CP_X_MAX, par.DZ_X_MIN, par.DZ_SEED_MIN, par.reseedSlope, par.reseedB,
        par.FRAP_t, par.FRAP_x1, par.FRAP_x2, par.FRAP_y1, par.FRAP_y2, par.no_grf, par.outputJPEG, par.blurRadiusJPEG);

    php.SetNumBins(par.bin_N);

    if (par.NoMature_mode)
    {
        php.SetNoMatureMode(true);
    }

    switch (init_type)
    {
    case INIT_TYPE_RND:
    {
        php.Init(par.frac_treadmill, par.frac_grow_out, par.frac_grow_in, par.theta_max, par.frac_seed_CP,
                 par.frac_seed_distal, par.frac_seed_middle, par.frac_seed_middle_slope, par.L0);
        break;
    }
    case INIT_TYPE_FIXED_L:
    {
        php.InitFixedLength(par.L0, par.theta_max);
        break;
    }
    default:
    {
        exit(0);
    }
    }

    if (par.no_grf == 0)
    {
        php.SetNoGRF();
    }

    php.RunSync(par.t_min, par.t_max, par.rnd_seed);

    return 0;
}
