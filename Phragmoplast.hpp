/**
 * @file Phragmoplast.hpp
 * @brief Phragmoplast Class Declaration
 * @version 2
 * @date 2021-10-25
 * 
 */

#ifndef PHP_SIM_HPP
#define PHP_SIM_HPP

#include "php_const.hpp"

#include "Bitmap.hpp"
#include "FrameName.hpp"
#include "Microtubule.hpp"
#include "RandomGen.hpp"



using namespace std;

class Microtubule;

//
// Phragmoplast made out of MT
//

/**
 * @brief Phragmoplast made out of microtubules; Inherits RandomGen and FrameName
 * 
 */
class Phragmoplast : public RandomGen, public FrameName
{
private:
    double time;                       // Current time
    double timeMax;                    // Max time
    int numPixelX;                     // Number pixels horizontaly
    int numPixelY;                     // Number pixels verticaly
    double figureStep;                 // time step between figures
    double MaxSupperposition_;         // Max supperposition in images
    double ratePolym;                  // micron/s
    double rateDepolym;                // shrinking rate in micron/s
    double rateGrowth2Shrink;          // catastrophe rate in 1/s
    double rateGrowth2ShrinkMature;    // catastrophe rate when fully grown in 1/s
    double rateGrowth2ShrinkMidzone;   // catastrophe rate when in the midzone in 1/s
    double ratePause2Growth;           // Rate pause -> growing
    double ratePause2Shrink;           // Rate pause -> shrinking
    double rateGrowth2Pause;           // Rate growing -> pause
    double rateShrink2Pause;           // Rate shrinking -> pause
    double ratePause2ShrinkCellPlate;  // Rate shrinking -> pause at CP
    double ratePause2ShrinkDistalZone; // Rate shrinking -> pause at DZ
    double rateShrink2Growth;          // recovery rate in 1/s

    // Rates for -end
    double r_me_polym_;   // Polymerisation rate
    double r_me_depolym_; // Depolymerisation rate
    double r_me_gs_;      // Catastrophe rate
    double r_me_sg_;      // Rescure rate
    double r_me_pg_;      // Rate pause -> growing
    double r_me_ps_;      // Rate pause -> shrinking
    double r_me_gp_;      // Rate growing -> pause
    double r_me_sp_;      // Rate shrinking -> pause

    double rateReseed;      // Reseeding rate in 1/s
    double reseedSlope;     // Slope of the reseed line
    double reseedIntercept; // Y-intercept of the reseed line

    int distributionType;  // Type of seeding distribution
    double seedParameterA; // Parameter A for the seeding distribution
    double seedParameterB; // Parameter B for the seeding distribution

    std::string figureName;     // Name of the figure
    std::string kymographName;  // Name of the kymograph  file
    std::string luminosityName; // Name of the luminosity file
    std::string histogramName;  // Name of the histrogram file

    Bitmap *kymographBitmap; // Kymograph
    int kymographIndex;      // Current number of kymograph lines
    int kymographLines;      // The number of kymograph lines
    double kymographAmp;     // Kymograph luminosity
    int numberBins;          // Number of bins in histogram

    vector<Microtubule> MTArray; // vector of MT

    double densityTubulin; // Length (um) density of MT per um^2

    int gloablSimIndex; // SimIndex for folder ID in large simulation runs

    double initialLengthMT; // Initial total length of MT
    double totalLengthMT;   // Total length of MT

public:
    double CellPlateXMax;  // Xmax boundary of the cell plate region
    double DistalZoneXMin; // Xmin boundary of the distal region
    double DistalSeedMin;  // Xmin boundary when seeding of the distal region
    double InitialLength;  // Maximum initial length or initial length (fixed init)
    double ThicknessMT;    // Thickness of phragmoplast (length of MT in microns)
    double WidthMT;        // Width of phragmoplast (radial width) (microns)
    double DeltaX;         // MT step growth (micron, use 100 nm?)
    double ThetaMax;       // Maximum theta value (degrees)
    int NumberMT;          // Number of microtubules

    //
    // FRAP parameters
    //

    double FRAPTime;   // Time FRAP occurs
    double FRAPX1;     // FRAP region x1
    double FRAPY1;     // FRAP region y1
    double FRAPX2;     // FRAP region x2
    double FRAPY2;     // FRAP region y2
    bool NoMatureMode; // No mature mode enable = true

    bool OutputBMA;        // If true: generate BMA files
    bool OutputJPEG;       // If true: generate JPEG files
    double JPEGBlurRadius; // Blur radius for JPEG images

    Phragmoplast();
    Phragmoplast(int N_MT, double Thickness, double Width, double dx, int Npixel_X, int Npixel_Y, double grf_dt,
                 double r_polym, double r_depolym, double r_gs, double r_sg, double r_reseed, double r_pg, double r_ps,
                 double r_gp, double r_sp, double r_ps_CP, double r_ps_DZ, double r_me_polym, double r_me_depolym,
                 double r_me_gs, double r_me_sg, double r_me_pg, double r_me_ps, double r_me_gp, double r_me_sp,
                 const char *figname, double grf_max, double kymograph_amp, double dens_MT, int gloablSimulationIndex, double t_max,
                 double seedPA, double seedPB, int disttype, double CPMax, double DZMin, double dzSeedMin, double slope,
                 double b, double frapT, double x1, double x2, double y1, double y2, bool outputBma, bool outputJpeg,
                 double jpegBlurRad);
    ~Phragmoplast();

    void Init(double frac_treadmill, double frac_grow_out, double frac_grow_in, double theta_max, double frac_seed_CP,
              double frac_seed_distal, double frac_seed_middle, double frac_seed_middle_slope, double l0);

    void InitFixedLength(double initialLength, double theta_max);

    void RunSync(double t_min, double t_max, long seed = -1);

    void SyncAllMTUntil(double tmax);

    void Bleach();

    double TotalLength(double &average_L);

    inline double Random(double xmax)
    {
        return (next_0_1() * xmax);
    };
    inline double RandomTime(double rate)
    {
        return (-log(Next()) / rate);
    };

    void MakeFigure(double t);
    void MakeHistogram(double t);
    void SaveList(double t);
    void AddKymographLine(Bitmap &bm);
    void MakeKymograph();

    inline void SetNoMatureMode(bool on)
    {
        NoMatureMode = on;
    }
    inline void SetNoGRF()
    {
        OutputBMA = false;
    };
    inline void SetNumBins(int bin_N)
    {
        numberBins = bin_N;
    };

    inline void IncrementTotalLength(double L)
    {
        totalLengthMT += L;
    };
    inline void DecrementTotalLength(double L)
    {
        totalLengthMT -= L;
    };
    inline double GetTotalLength(void)
    {
        return totalLengthMT;
    };
    inline double RelativeMT()
    {
        return (double)(totalLengthMT / initialLengthMT);
    };

    inline double GetThetaMax()
    {
        return (ThetaMax);
    };
    inline double GetRateDepolym()
    {
        return (rateDepolym);
    };

    inline double GetRatePause2Growth()
    {
        return (ratePause2Growth);
    };
    inline double GetRatePause2Shrink()
    {
        return (ratePause2Shrink);
    };

    inline double GetRateShrink2Pause()
    {
        return (rateShrink2Pause);
    };
    inline double GetRateShrink2PauseCellPlate()
    {
        return (ratePause2ShrinkCellPlate);
    };
    inline double GetRateShrink2PauseDistalZone()
    {
        return (ratePause2ShrinkDistalZone);
    };

    inline double get_r_me_polym()
    {
        return (r_me_polym_);
    };
    inline double get_r_me_depolym()
    {
        return (r_me_depolym_);
    };
    inline double get_r_me_gs()
    {
        return (r_me_gs_);
    };
    inline double get_r_me_sg()
    {
        return (r_me_sg_);
    };
    inline double get_r_me_pg()
    {
        return (r_me_pg_);
    };
    inline double get_r_me_ps()
    {
        return (r_me_ps_);
    };
    inline double get_r_me_gp()
    {
        return (r_me_gp_);
    };
    inline double get_r_me_sp()
    {
        return (r_me_sp_);
    };

    inline double RateShrink2Growth()
    {
        return (rateShrink2Growth * RelativeMT());
    };
    inline double RateShrink2Pause()
    {
        return (rateShrink2Pause * RelativeMT());
    };
    inline double RatePause2Growth()
    {
        return (ratePause2Growth * RelativeMT());
    };
    inline double RatePause2Shrink()
    {
        return (ratePause2Shrink / (RelativeMT() + 1e-6));
    };
    inline double RatePause2ShrinkCellPlate()
    {
        return (ratePause2ShrinkCellPlate / (RelativeMT() + 1e-6));
    };
    inline double RatePause2ShrinkDistalZone()
    {
        return (ratePause2ShrinkDistalZone / (RelativeMT() + 1e-6));
    };

    inline double R_me_polym()
    {
        return (r_me_polym_ * RelativeMT());
    };
    inline double R_me_sg()
    {
        return (r_me_sg_ * RelativeMT());
    };
    inline double R_me_pg()
    {
        return (r_me_pg_ * RelativeMT());
    };
    inline double R_me_ps()
    {
        return (r_me_ps_ / (RelativeMT()));
    };
    inline double R_me_sp()
    {
        return (r_me_sp_ * RelativeMT());
    };

    // MT dependant rates
    inline double RateGrowth2Pause(void)
    {
        return rateGrowth2Pause / (RelativeMT() + 1e-6);
    };
    double RateGrowth2Shrink(void)
    {
        return rateGrowth2Shrink / (RelativeMT() + 1e-6);
    };
    inline double R_me_gs(void)
    {
        return r_me_gs_ / (RelativeMT() + 1e-6);
    };
    inline double R_me_gp(void)
    {
        return r_me_gp_ / (RelativeMT() + 1e-6);
    };
    inline double RatePolym()
    {
        return (ratePolym * RelativeMT());
    };
    double RateReseed(double xInit, int type);
};
#endif
