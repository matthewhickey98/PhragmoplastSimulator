#include "Phragmoplast.hpp"

using namespace std;


/**
 * @brief Construct a new Phragmoplast
 *
 * @param N_MT                      number of microtubles
 * @param Thickness                 thickness in microns
 * @param Width                     width in microns
 * @param dx                        step in x change in microns
 * @param Npixel_X                  number of pixels in X
 * @param Npixel_Y                  number of pixels in Y
 * @param grf_dt                    time step in seconds between images
 * @param r_polym                   rate of polymerization
 * @param r_depolym                 rate of depolymerization
 * @param r_gs                      rate of growth to shrink
 * @param r_sg                      rate of shrink to grow
 * @param r_reseed                  rate of reseeding
 * @param r_pg                      rate of pause to grow
 * @param r_ps                      rate of pause to shrink
 * @param r_gp                      rate of grow to pause
 * @param r_sp                      rate of shrink to pause
 * @param r_ps_CP                   rate of pause to shrink in cell plate region
 * @param r_ps_DZ                   rate of pause to shrink in distal zone region
 * @param r_me_polym                rate of - end polymerization
 * @param r_me_depolym              rate of - end depolymerization
 * @param r_me_gs                   rate of - end growth to shrink
 * @param r_me_sg                   rate of - end shrink to grow
 * @param r_me_pg                   rate of - end pause to grow
 * @param r_me_ps                   rate of - end pause to shrink
 * @param r_me_gp                   rate of - end grow to pause
 * @param r_me_sp                   rate of - end shrink to pause
 * @param figname                   name of the output file
 * @param grf_max                   maximum number of images
 * @param kymograph_amp             amplitude of the kymograph
 * @param dens_MT                   density of tublulin
 * @param gloablSimulationIndex     index of the simulation
 * @param t_max                     maximum time in seconds
 * @param seedPA                    seed parameter A in the distribution
 * @param seedPB                    seed parameter B in the distribution
 * @param disttype                  type of distribution
 * @param CPMax                     maximum X coordinate of the cell plate region
 * @param DZMin                     minimum X coordinate of the distal zone region
 * @param dzSeedMin                 mimum X coordinate of the distal zone seeding
 * @param Slope                     slope of the distal zone seeding
 * @param B                         intercept of the distal zone seeding
 * @param frapT                     time of the frap in seconds`
 * @param x1                        x coordinate of the first point of the frap
 * @param x2                        x coordinate of the second point of the frap
 * @param y1                        y coordinate of the first point of the frap
 * @param y2                        y coordinate of the second point of the frap
 * @param outputBma                 output BMA files
 * @param outputJpeg                output JPEG files
 * @param jpegBlurRad               blur radius of the JPEG files
 */
Phragmoplast::Phragmoplast(int N_MT, double Thickness, double Width, double dx, int Npixel_X, int Npixel_Y,
                           double grf_dt, double r_polym, double r_depolym, double r_gs, double r_sg, double r_reseed,
                           double r_pg, double r_ps, double r_gp, double r_sp, double r_ps_CP, double r_ps_DZ,
                           double r_me_polym, double r_me_depolym, double r_me_gs, double r_me_sg, double r_me_pg,
                           double r_me_ps, double r_me_gp, double r_me_sp, const char *figname, double grf_max,
                           double kymograph_amp, double dens_MT, int gloablSimulationIndex, double t_max, double seedPA, double seedPB,
                           int disttype, double CPMax, double DZMin, double dzSeedMin, double Slope, double B,
                           double frapT, double x1, double x2, double y1, double y2, bool outputBma, bool outputJpeg,
                           double jpegBlurRad)
    : NumberMT(N_MT), ThicknessMT(Thickness), WidthMT(Width), time(0.0), DeltaX(dx), numPixelX(Npixel_X),
      numPixelY(Npixel_Y), figureStep(grf_dt), MaxSupperposition_(grf_max), ratePolym(r_polym), rateDepolym(r_depolym),
      rateGrowth2Shrink(r_gs), rateShrink2Growth(r_sg), rateReseed(r_reseed), ratePause2Growth(r_pg),
      ratePause2Shrink(r_ps), rateGrowth2Pause(r_gp), rateShrink2Pause(r_sp), ratePause2ShrinkCellPlate(r_ps_CP),
      ratePause2ShrinkDistalZone(r_ps_DZ), r_me_polym_(r_me_polym), r_me_depolym_(r_me_depolym), r_me_gs_(r_me_gs),
      r_me_sg_(r_me_sg), r_me_pg_(r_me_pg), r_me_ps_(r_me_ps), r_me_gp_(r_me_gp),
      r_me_sp_(r_me_sp), FrameName::FrameName(figname, ".bma"), kymographBitmap(0), kymographAmp(kymograph_amp),
      OutputBMA(outputBma), numberBins(100), densityTubulin(dens_MT), gloablSimIndex(gloablSimulationIndex), seedParameterA(seedPA),
      seedParameterB(seedPB), distributionType(disttype), CellPlateXMax(CPMax), DistalZoneXMin(DZMin),
      DistalSeedMin(dzSeedMin), reseedSlope(Slope), reseedIntercept(B), FRAPTime(frapT), FRAPX1(x1), FRAPX2(x2),
      FRAPY1(y1), FRAPY2(y2), OutputJPEG(outputJpeg), JPEGBlurRadius(jpegBlurRad)
{
    NoMatureMode = false;
    initialLengthMT = totalLengthMT = densityTubulin * ThicknessMT * WidthMT;
    figureName.resize(1024);
    stringstream ss;
    ss << figname << "_N" << NumberMT << "_dMT" << densityTubulin << "_X" << Thickness << "_Y" << Width << "_dx" << dx
       << "_r+" << r_polym << "_r-" << r_depolym << "_Rgs" << r_gs << "_Rsg" << r_sg << "_Rreseed" << r_reseed << "_Rpg"
       << r_pg << "_Rps" << r_ps << "_Rgp" << r_gp << "_Rsp" << r_sp << "_RspCP" << r_ps_CP << "_RspDZ" << r_ps_DZ
       << "_";
    figureName = ss.str();
    SetPrefix(figureName);
    cout << "Phragmoplast figureName=" << figureName << "\n";
    MTArray.resize(NumberMT);
    kymographLines = (int)floor(t_max / figureStep + 0.999) + 1;
    kymographBitmap = new Bitmap(Npixel_X, kymographLines, Thickness, Width, t_max, 1.0);
    cout << "Slope:" << reseedSlope << "    b: " << reseedIntercept << endl;
}
Phragmoplast::~Phragmoplast()
{
    delete kymographBitmap;
}

/**
 * @brief Initialize a new simulation with random length
 *
 * @param frac_treadmill            fraction of the MTs that treadmill
 * @param frac_grow_out             fraction of the MTs that grow out
 * @param frac_grow_in              fraction of the MTs that grow in
 * @param theta_max                 maximum angle of the MTs
 * @param frac_seed_CP              fraction of the MTs that seed in the cell plate
 * @param frac_seed_distal          fraction of the MTs that seed in the distal zone
 * @param frac_seed_middle          fraction of the MTs that seed in the middle
 * @param frac_seed_middle_slope    slope of the middle seed
 * @param l0                        initial length of the MTs
 */
void Phragmoplast::Init(double frac_treadmill, double frac_grow_out, double frac_grow_in, double theta_max,
                        double frac_seed_CP, double frac_seed_distal, double frac_seed_middle,
                        double frac_seed_middle_slope, double l0)
{
    double l_treadmill, l_grow_out, l_grow_in, tot;
    double l_seed_CP, l_seed_distal, l_seed_middle, p;
    double x1, y1, L, theta;

    int type;
    InitialLength = l0;
    ThetaMax = theta_max;
    tot = frac_treadmill + frac_grow_out + frac_grow_in;
    frac_treadmill /= tot;
    frac_grow_out /= tot;
    frac_grow_in /= tot;
    l_treadmill = frac_treadmill;
    l_grow_out = frac_grow_out + l_treadmill;
    l_grow_in = frac_grow_in + l_grow_out;
    // tot = frac_seed_CP + frac_seed_distal + frac_seed_middle;
    // frac_seed_CP /= tot;
    // frac_seed_distal /= tot;
    // frac_seed_middle /= tot;
    l_seed_CP = frac_seed_CP;
    l_seed_distal = frac_seed_distal + l_seed_CP;
    l_seed_middle = frac_seed_middle + l_seed_distal;

    std::ofstream initFile(figureName + "_InitFile.csv");
    initFile << "Index,x1,y1,L,theta,type,location" << std::endl;

    cout << frac_seed_middle * NumberMT << endl;

    int numDistal = 0;

    std::random_device rd;
    std::cauchy_distribution<double> cauchyDist(seedParameterA, seedParameterB);
    std::normal_distribution<double> normDist(seedParameterA, seedParameterB);
    std::weibull_distribution<double> weibullDist(seedParameterA, seedParameterB);
    std::exponential_distribution<double> exponentialDist(seedParameterA);

    for (int i = 0; i < NumberMT; ++i)
    {
        L = Random(InitialLength); // NEW thicknessMT -> initialLength
        y1 = Random(WidthMT);
        theta = Random(theta_max);
        p = Random(1);
        if (p < l_treadmill)
        {
            type = PHP_TYPE_TREADMILL;
        }
        else if (p < l_grow_out)
        {
            type = PHP_TYPE_GROW_OUTSIDE;
        }
        else
        {
            type = PHP_TYPE_GROW_INSIDE;
        }

        //
        // cellplate and middle will be seeded on gradient distribution. The rest
        // of the seeds will be randomly seeded in the distal region
        //

        if (numDistal < (int)(frac_seed_distal * NumberMT))
        {
            x1 = (ThicknessMT - DistalSeedMin) * Random(1) + DistalSeedMin;

            type = PHP_TYPE_GROW_INSIDE;
            initFile << i << "," << x1 << "," << y1 << "," << L << "," << theta << "," << type << ",D" << std::endl;
            numDistal++;
            cout << "NumDistal: " << numDistal << endl;
        }
        else
        {
            do
            {
                if (distributionType == 0)
                {
                    x1 = cauchyDist(rd);
                }

                else if (distributionType == 1)
                {
                    x1 = normDist(rd);
                }

                else if (distributionType == 2)
                {
                    x1 = weibullDist(rd);
                }

                else
                {
                    x1 = DistalSeedMin * (1 - exponentialDist(rd));
                }

            } while (x1 < 0 || x1 > DistalSeedMin);

            if (x1 < CellPlateXMax)
            {
                type = PHP_TYPE_GROW_OUTSIDE;
            }

            if (x1 > DistalZoneXMin)
            {
                type = PHP_TYPE_GROW_INSIDE;
            }

            initFile << i << "," << x1 << "," << y1 << "," << L << "," << theta << "," << type << ",M" << std::endl;
        }
        MTArray[i].Init(x1, y1, L, theta, type, this);
    }
    initFile.close();
}

/**
 * @brief Initialize a new simulation with fixed length
 *
 * @param initialLength
 * @param theta_max
 */
void Phragmoplast::InitFixedLength(double initialLength, double theta_max)
{
    double y1, theta;
    std::cout << "init_fixed_L(" << initialLength << "," << theta_max << ")\n";
    ThetaMax = theta_max;
    for (int i = 0; i < NumberMT; ++i)
    {
        y1 = this->Next() * WidthMT;
        theta = (this->Next() - 0.5) * 2.0 * theta_max;
        MTArray[i].Init(ThicknessMT, y1, initialLength, theta, PHP_TYPE_GROW_INSIDE, this);
        MTArray[i].SetPlusState(PHP_STATE_SHRINKING);
    }
}

/**
 * @brief Evolve all MT synchonously from t_min to t_max
 *
 * @param t_min     start time
 * @param t_max     end time
 * @param seed      seed for the random number generator
 */
void Phragmoplast::RunSync(double t_min, double t_max, long seed)
{
    timeMax = t_max;
    double next_grf_t;
    ofstream ofs;
    bool just_frapped = false;
    stringstream ss;
    ss << "DATA_" << this->gloablSimIndex;
    string folderName;
    ss >> folderName;
    std::cout << folderName << endl;
    ofstream lengthsFile(figureName + "_LengthsFile.csv");
    ofstream fractionFile(figureName + "_Fractions.csv");
    ofstream relMTFile(figureName + "_RelativeMTDensity.csv");
    ofstream totalLengthFile(figureName + "_TotalMTLength.csv");

    std::cout << "run_sync(" << t_min << "," << t_max << "," << seed << ")\n";

    kymographIndex = 0;
    kymographLines = (int)floor(t_max / figureStep + 0.999) + 1;
    delete kymographBitmap;
    std::cout << "FRAPx1__=" << FRAPX1 << " FRAPy1__=" << FRAPY1 << "FRAPx2__=" << FRAPX2 << " FRAPy2__=" << FRAPY2
              << "\n";

    kymographBitmap = new Bitmap(numPixelX, kymographLines, ThicknessMT, WidthMT, t_max, 1.0);
    if (kymographAmp <= 0)
        kymographBitmap->set_Max(-1);

    kymographName = figureName + "_kymograph.jpeg";

    luminosityName = figureName + "_luminosity.txt";
    histogramName = figureName + "_histogram";

    ofs.open(luminosityName);
    ofs << "t,luminosity,numEmpty,numMature,numGrowing,numShrinking,numPause,growOut,growIn,treadmill,totalLen,averageLen,relMT\n";
    ofs.close();

    for (int i = 0; i < NumberMT; ++i)
    {
        MTArray[i].SetTime(t_min);
    }

    if ((time = t_min) < 1e-16)
    {
        SyncAllMTUntil(0.0);
        time = 0.0;
    }

    MakeFigure(time);
    MakeHistogram(time);
    next_grf_t = time + figureStep;

    fractionFile
        << "Time, SeededCount, UnseededCount, [0, 0.06], (0.06, 0.12], (0.12, 0.18], (0.18, 0.24], (0.24, 0.36], (0.36, 0.42], \
     (0.42, 0.48], (0.48, 0.54], (0.54, 0.60], (0.66, 0.72], (0.72, 0.78], (0.78, 0.84], (0.84, 0.90], (0.96, 1.02], (1.02, 1.08], \
     (1.08, 1.14], (1.14, 1.20], (1.20, 1.26], (1.26, 1.32], (1.32, 1.38], (1.38, 1.42], (1.42, 1.5]"
        << endl;

    while (time <= t_max)
    {
        //
        // Calculate the fractions and write them to a file
        //

        double averageL;
        totalLengthFile << time << "," << TotalLength(averageL) << "," << averageL << endl;
        lengthsFile << time << ",";
        relMTFile << time << "," << RelativeMT() << "," << totalLengthMT << "," << initialLengthMT << endl;

        int count[25] = {0};
        int emptyCount = 0;
        int seededCount = 0;
        for (int i = 0; i < NumberMT; i++)
        {
            double length = MTArray[i].Length();
            lengthsFile << length << ",";
            if (length <= 0)
            {
                emptyCount++;
            }
            else
            {
                seededCount++;

                for (int j = 0; j < 25; j++)
                {
                    if (length <= 0.06 + 0.06 * j)
                    {
                        count[j]++;
                        break;
                    }
                }
            }
        }
        lengthsFile << endl;

        fractionFile << time << ", " << seededCount << ", " << emptyCount << ", ";
        for (int i = 0; i < 25; i++)
        {
            fractionFile << count[i] << ", ";
        }

        fractionFile << endl;

        if ((FRAPTime > 0) && (FRAPTime <= next_grf_t))
        {
            SyncAllMTUntil(FRAPTime);
            Bleach();
            std::cout << " FRAP at t=" << time << "\n";
            FRAPTime = -1;
            just_frapped = true;
        }

        if (just_frapped)
        {
            just_frapped = false;
            if (time <= next_grf_t - figureStep * 0.01)
            {
                SyncAllMTUntil(next_grf_t);
            }
        }
        else
        {
            SyncAllMTUntil(next_grf_t);
        }

        MakeFigure(next_grf_t);
        MakeHistogram(next_grf_t);

        std::cout << "t = " << time << "\n";
        next_grf_t += figureStep;
    }
    if (time > next_grf_t - figureStep * 0.5)
    {
        MakeFigure(next_grf_t);
        MakeHistogram(next_grf_t);
    };

    fractionFile.close();
    lengthsFile.close();
    relMTFile.close();
    totalLengthFile.close();
    MakeKymograph();

    std::cout << "Run finished: t = " << time << "\n";
}

/**
 * @brief Syn all microtubules until a certain time
 *
 * @param tmax  the time to sync to
 */
void Phragmoplast::SyncAllMTUntil(double tmax)
{
    double dt, dt_min;
    int MT_min = -1, event, event_min;

    while (time < tmax)
    {
        dt_min = 1e64;

        for (int i = 0; i < NumberMT; ++i)
        {
            dt = MTArray[i].SelectNextEvent();

            if (dt < dt_min)
            {
                dt_min = dt;
                MT_min = i;
                event_min = MTArray[i].GetMinEvent();
            }
        }
        if (time + dt_min > tmax) // time to stop
        {
            time = tmax;
            return;
        }
        MTArray[MT_min].DoEvent(event_min);
        time += dt_min;

#ifdef MONITOR
        if (trace++ > 10)
        {
            cout << "t=" << t_ << " last event:" << event_min << "\n";
            trace = 0;
            for (int i = PHP_EVENT_P_GROW; i <= PHP_EVENT_P_RESEED; i++) // MONITOR
            {
                cout << "event_no[" << i << "]:" << event_no[i] << "\n";
                event_no[i] = 0;
            }
            for (int i = PHP_EVENT_M_GROW; i <= PHP_EVENT_M_RESEED; i++) // MONITOR
            {
                cout << "event_no[" << i << "]:" << event_no[i] << "\n";
                event_no[i] = 0;
            }
        }
#endif
        // cout << "t_ = " << t_ << "    tmax = " << tmax << endl;
    }
};

/**
 * @brief Bleach each microtubule
 *
 */
void Phragmoplast::Bleach()
{
    for (int i = 0; i < NumberMT; ++i)
    {
        MTArray[i].Bleach();
    }
}

//
// Compute the total MT length and average length
// Write average length in double &average_L
//

/**
 * @brief Compute the total and average length of all the microtubules
 *
 * @param average_L     the average length of all the microtubules
 * @return double       the total length of all the microtubules
 */
double Phragmoplast::TotalLength(double &average_L)
{
    double L = 0.0;
    int n = 0;

    for (int i = 0; i < NumberMT; ++i)
    {
        if (MTArray[i].GetPlusState() != PHP_STATE_EMPTY)
        {
            ++n;
            L += MTArray[i].Length();
        }
    }
    average_L = (n > 0) ? L / n : 0;
    return (L);
}

void Phragmoplast::MakeFigure(double t)
{
    Bitmap bm(numPixelX, numPixelY, ThicknessMT, WidthMT, t, MaxSupperposition_);
    ofstream ofs;
    double empty = 0, mature = 0, growing = 0, shrinking = 0, pause = 0, L, average_L;
    double Grow_out = 0, Grow_in = 0, teadmill = 0;

    for (int i = 0; i < NumberMT; ++i)
    {
        MTArray[i].Draw(bm);

        if (MTArray[i].GetPlusState() == PHP_STATE_EMPTY)
        {
            ++empty;
        }
        if (MTArray[i].GetPlusState() == PHP_STATE_MATURE)
        {
            ++mature;
        }
        if (MTArray[i].GetPlusState() == PHP_STATE_GROWING)
        {
            ++growing;
        }
        if (MTArray[i].GetPlusState() == PHP_STATE_SHRINKING)
        {
            ++shrinking;
        }
        if (MTArray[i].GetPlusState() == PHP_STATE_PAUSE)
        {
            ++pause;
        }
        if (MTArray[i].GetType() == PHP_TYPE_GROW_OUTSIDE)
        {
            ++Grow_out;
        }
        if (MTArray[i].GetType() == PHP_TYPE_GROW_INSIDE)
        {
            ++Grow_in;
        }
        if (MTArray[i].GetType() == PHP_TYPE_TREADMILL)
        {
            ++teadmill;
        }
    }
    // bm.gaussian_blur(blur_radius_);
    AddKymographLine(bm);
    // bm.output(fname_.next_name());
    NextName();
    if (OutputBMA == 1)
    {
        bm.Save(CurrentName());
        SaveList(t);
    }
    if (OutputJPEG == 1)
    {
        if (JPEGBlurRadius > 0.2)
        {
            bm.GaussianBlur(JPEGBlurRadius);
        }
        bm.Output(CurrentNameOtherSuffix(".jpg"), true, 0);
    }

    L = TotalLength(average_L);
    ofs.open(luminosityName, std::ios::app);
    ofs << t << "," << bm.Luminosity(FRAPX1 + DeltaX, FRAPY1 + DeltaX, FRAPX2 - DeltaX, FRAPY2 - DeltaX) / NumberMT
        << "," << empty << "," << mature << "," << growing << "," << shrinking << "," << pause << "," << Grow_out << ","
        << Grow_in << "," << teadmill << "," << L << "," << average_L << "," << RelativeMT() << "\n";
    ofs.close();
}

/**
 * @brief Make a histogram of the MT lengths
 *
 * @param t     the current time
 */
void Phragmoplast::MakeHistogram(double t)
{
    ofstream ofs;
    double hist[1000];
    int j, N;
    char filename[PHP_PATH_LENGTH + 100];
    // cout<< "make_histogram("<<t<<")\n";
    sprintf(filename, "%s_t%g.txt", histogramName.c_str(), t);
    for (j = 0; j < numberBins; ++j)
    {
        hist[j] = 0;
    }
    N = 0;

    for (int i = 0; i < NumberMT; ++i)
    {
        if (MTArray[i].GetPlusState() != PHP_STATE_EMPTY)
        {
            j = (int)floor(MTArray[i].Length() * (numberBins - 1) / ThicknessMT);
            if (j >= numberBins)
                j = numberBins - 1;
            hist[j] += 1;
            ++N;
        }
    }

    ofs.open(filename, std::ios::out);
    ofs << "# t=" << t << "\n";
    ofs << "# N=" << N << "\n";
    for (j = 0; j < numberBins; ++j)
    {
        ofs << j * ThicknessMT / (numberBins - 1.0) << " " << hist[j] / N << "\n";
    }

    ofs.close();
}

/**
 * @brief Save a list of information about microtubule positions
 *
 * @param t         the current time
 */
void Phragmoplast::SaveList(double t)
{
    ofstream ofs(this->CurrentNameOtherSuffix(".txt"), std::ios::app);

    ofs << "#index type state x1 y1 x2 y2 bleached_x1 bleached_x2 \n";
    ofs << "#t=" << t << "\n";
    for (int i = 0; i < NumberMT; ++i)
    {
        ofs << i << " ";
        MTArray[i].TextSave(ofs);
    }
}

/**
 * @brief Add a line to the kymograph
 *
 * @param bm    the bitmap
 */
void Phragmoplast::AddKymographLine(Bitmap &bm)
{
    std::vector<double> v(numPixelX);

    bm.Average(v, kymographAmp);

    kymographBitmap->AddLine(v, kymographIndex++);
}

/**
 * @brief Make a kymograph
 *
 */
void Phragmoplast::MakeKymograph()
{
    kymographBitmap->Save(kymographName);
    kymographBitmap->Output(kymographName, false, this->time);
}

/**
 * @brief Get the rate of reseed, which is influenced by the position and type of the microtubule
 *
 * @param xInit     the seed x position
 * @param type      the type of the microtubule
 * @return double   the rate of reseed
 */
double Phragmoplast::RateReseed(double xInit, int type)
{
    double r_reseed = (rateReseed * (RelativeMT()));
    if (type == PHP_TYPE_GROW_INSIDE)
        r_reseed = (xInit * reseedSlope + reseedIntercept) * (RelativeMT());
    if (r_reseed > rateReseed)
        return rateReseed * (RelativeMT());
    else if (r_reseed < 0)
        return 0;
    else
        return r_reseed;
}
