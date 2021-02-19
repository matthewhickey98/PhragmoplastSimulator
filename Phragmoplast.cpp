#include "Phragmoplast.h"
using namespace std;

//************************************************************/
// Phragmoplast make out of MT
Phragmoplast::Phragmoplast(int N_MT, double Thickness, double Width, double dx,
                           int Npixel_X, int Npixel_Y, double grf_dt,
                           double r_polym, double r_depolym,
                           double r_gs,
                           double r_sg, double r_reseed,
                           double r_pg, double r_ps, double r_gp, double r_sp,
                           double r_ps_CP, double r_ps_DZ,
                           double r_me_polym, double r_me_depolym, double r_me_gs,
                           double r_me_sg,
                           double r_me_pg, double r_me_ps, double r_me_gp, double r_me_sp,
                           const char *figname,
                           double grf_max, double kymograph_amp, double dens_MT,
                           double macet_dz_len, double macet_mz, double macet_dz,
                           double macet_MT_len) : N_MT_(N_MT), Thickness_(Thickness), Width_(Width), dx_(dx), t_(0.0),
                                                  Npixel_X_(Npixel_X), Npixel_Y_(Npixel_Y), grf_dt_(grf_dt),
                                                  MaxSupperposition_(grf_max),
                                                  r_polym_(r_polym), r_depolym_(r_depolym),
                                                  r_gs_(r_gs), r_sg_(r_sg), r_reseed_(r_reseed),
                                                  r_pg_(r_pg), r_ps_(r_ps), r_gp_(r_gp), r_sp_(r_sp), r_ps_CP_(r_ps_CP),
                                                  r_ps_DZ_(r_ps_DZ),
                                                  r_me_polym_(r_me_polym), r_me_depolym_(r_me_depolym), r_me_gs_(r_me_gs),
                                                  r_me_sg_(r_me_sg),
                                                  r_me_pg_(r_me_pg), r_me_ps_(r_me_ps), r_me_gp_(r_me_gp), r_me_sp_(r_me_sp),
                                                  FRAP_t_(-1),
                                                  fname_(figname, ".bma"), kymograph_bm_(0), kymograph_amp_(kymograph_amp),
                                                  no_grf_(false), bin_N_(100), dens_MT_(dens_MT), macet_dz_len_(macet_dz_len),
                                                  macet_mz_(macet_mz), macet_dz_(macet_dz), macet_MT_len_(macet_MT_len)
{
    Microtuble::set_bounds(Thickness, Width, dx);
    tot_MT_0_ = tot_MT_ = dens_MT_ * Thickness_ * Width_;

    sprintf(fig_name_, "%s_N%d_dMT%g_X%g_Y%g_dx%g_r+%g_r-%g_rc%g_rR%g_rs%g_pg%g_ps%g_gp%g_sp%g_spCP%g_spDZ%g",
            figname,
            N_MT_, dens_MT_, Thickness, Width, dx, r_polym, r_depolym, r_gs,
            r_sg, r_reseed, r_pg, r_ps, r_gp, r_sp, r_ps_CP, r_ps_DZ);
    fname_.set_prefix(fig_name_);
    //cerr<<"Phragmoplast fig_name_="<<fig_name_<<"\n";
    mt_array_ = new Microtuble[N_MT_];
}

void Phragmoplast::init(double frac_treadmill, double frac_grow_out,
                        double frac_grow_in, double theta_max,
                        double frac_seed_CP, double frac_seed_distal,
                        double frac_seed_middle, double frac_seed_middle_slope, double L0)
{
    double l_treadmill, l_grow_out, l_grow_in, tot;
    double l_seed_CP, l_seed_distal, l_seed_middle, p;
    double x1, y1, L, theta;

    int type;

    theta_max_ = theta_max;
    tot = frac_treadmill + frac_grow_out + frac_grow_in;
    frac_treadmill /= tot;
    frac_grow_out /= tot;
    frac_grow_in /= tot;
    l_treadmill = frac_treadmill;
    l_grow_out = frac_grow_out + l_treadmill;
    l_grow_in = frac_grow_in + l_grow_out;

    tot = frac_seed_CP + frac_seed_distal + frac_seed_middle;
    frac_seed_CP /= tot;
    frac_seed_distal /= tot;
    frac_seed_middle /= tot;
    l_seed_CP = frac_seed_CP;
    l_seed_distal = frac_seed_distal + l_seed_CP;
    l_seed_middle = frac_seed_middle + l_seed_distal;

    for (int i = 0; i < N_MT_; ++i)
    {
        p = rnd_.next();
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

        p = rnd_.next();
        if ((p < l_seed_CP) || (type == PHP_TYPE_TREADMILL))
        {
            x1 = 0;
            if (type == PHP_TYPE_GROW_INSIDE)
            {
                type = PHP_TYPE_GROW_OUTSIDE;
            }
        }
        else if (p < l_seed_distal)
        {
            x1 = Thickness_ - Microtuble::min_L__;
            type = PHP_TYPE_GROW_INSIDE;
        }
        else
        {
            if (frac_seed_middle_slope > 1e-6)
            {
                x1 = rnd_.trapezium(frac_seed_middle_slope) * Thickness_;
            }
            else
            {
                x1 = rnd_.next() * Thickness_;
            }
        }

        L = rnd_.next()*L0; // NEW Thickness_ -> L0
        tot_MT_ -= L;
        cout<<"totalMT0: "<<tot_MT_0_<<" total MT: "<<tot_MT_<<" L: "<<L<<endl;
        y1 = rnd_.next() * Width_;
        theta = (rnd_.next() - 0.5) * 2.0 * theta_max;

        mt_array_[i].init(x1, y1, L, theta, type, this);
    }
}

void Phragmoplast::init_fixed_L(double L0, double theta_max)
{
    double y1, theta;

    cerr << "init_fixed_L(" << L0 << "," << theta_max << ")\n";

    theta_max_ = theta_max;

    for (int i = 0; i < N_MT_; ++i)
    {
        y1 = rnd_.next() * Width_;
        theta = (rnd_.next() - 0.5) * 2.0 * theta_max;

        mt_array_[i].init(Thickness_, y1, L0, theta, PHP_TYPE_GROW_INSIDE, this);
        mt_array_[i].set_state_p(PHP_STATE_SHRINKING);
    }
}

// Evolve all MT synchonously from t_min to t_max
void Phragmoplast::run_sync(double t_min, double t_max, long seed)
{
    double next_grf_t;
    ofstream ofs;
    bool just_frapped = false;

    cerr << "run_sync(" << t_min << "," << t_max << "," << seed << ")\n";

    index_kymograph_ = 0;
    kymograph_lines_ = (int)floor(t_max / grf_dt_ + 0.999) + 1;
    delete kymograph_bm_;
    //cerr << "kymograph_lines_="<<kymograph_lines_<<"\n";
    //cerr << "Npixel_X_="<<Npixel_X_<<"\n";
    //cerr << "Thickness_="<<Width_<<"\n";
    cerr << "FRAPx1__=" << Microtuble::FRAPx1__ << " FRAPy1__=" << Microtuble::FRAPy1__ << "FRAPx2__=" << Microtuble::FRAPx2__ << " FRAPy2__=" << Microtuble::FRAPy2__ << "\n";

    kymograph_bm_ = new Bitmap(Npixel_X_, kymograph_lines_, Thickness_,
                               Width_, t_max, 1.0);
    if (kymograph_amp_ <= 0)
    {
        kymograph_bm_->set_Max(-1);
    }

    sprintf(kymograph_name_, "%s_kymograph.jpeg", fig_name_);
    sprintf(luminosity_name_, "%s_luminosity.txt", fig_name_);
    sprintf(histogram_name_, "%s_histogram", fig_name_);

    ofs.open(luminosity_name_);
    ofs << "#t luminosity empty mature growing shrinking pause gro gri tm tota_L average_L\n";
    ofs.close();
    for (int i = 0; i < N_MT_; ++i)
    {
        mt_array_[i].set_t(t_min);
    }

    cerr << "t_min = " << t_min << "\n";

    if ((t_ = t_min) < 1e-16)
    {
        evolve_sync_all_MT_until(0.0);
        t_ = 0.0;
    }
    ofstream lengthsFile("LengthsFile.csv");
    lengthsFile << endl << "T = " << t_ << endl;

    for (int i = 0; i < N_MT_; i++)
    {
        lengthsFile << mt_array_[i].L() << ","<<endl;
    }       
            
    make_fig(t_);
    make_histogram(t_);

    next_grf_t = t_ + grf_dt_;
    ofstream fractionFile("Fractions.csv");
    fractionFile << "Time, SeededCount, UnseededCount, 0_.1, .1_.2, .3_.4, .5_.6, .7_.8, .9_1.0, 1.1_1.2, 1.3_1.4, 1.4_1.5, 1.5+" << endl;
    
    while (t_ <= t_max)
    {
        //get the fractions
        //write fractions to Fractions.csv
        if ((int)t_ % 60 == 0 || (int)t_ == 0)
        {
            lengthsFile << endl
                        << "T = " << t_ << endl;

            int count[16] = {0};
            int emptyCount = 0;
            int seededCount = 0;
            for (int i = 0; i < N_MT_; i++)
            {
                double length = mt_array_[i].L();
                if (length <= 0)
                {
                    emptyCount++;
                }
                else
                {
                    seededCount++;
                    lengthsFile << length << endl;
                    if (length <= .1)
                    {
                        count[0]++;
                        continue;
                    }
                    else if (length <= .2)
                    {
                        count[1]++;
                        continue;
                    }
                    else if (length <= .3)
                    {
                        count[2]++;
                        continue;
                    }
                    else if (length <= .4)
                    {
                        count[3]++;
                        continue;
                    }
                    else if (length <= .5)
                    {
                        count[4]++;
                        continue;
                    }
                    else if (length <= .6)
                    {
                        count[5]++;
                        continue;
                    }
                    else if (length <= .7)
                    {
                        count[6]++;
                        continue;
                    }
                    else if (length <= .8)
                    {
                        count[7]++;
                        continue;
                    }
                    else if (length <= .9)
                    {
                        count[8]++;
                        continue;
                    }
                    else if (length <= 1.0)
                    {
                        count[9]++;
                        continue;
                    }
                    else if (length <= 1.1)
                    {
                        count[10]++;
                        continue;
                    }
                    else if (length <= 1.2)
                    {
                        count[11]++;
                        continue;
                    }
                    else if (length <= 1.3)
                    {
                        count[12]++;
                        continue;
                    }
                    else if (length <= 1.4)
                    {
                        count[13]++;
                        continue;
                    }
                    else if (length <= 1.5)
                    {
                        count[14]++;
                        continue;
                    }
                    else
                    {
                        count[15]++;
                        continue;
                    }
                }
            }

            double fracs[16];
            for (int i = 0; i < 16; i++)
                fracs[i] = (double)count[i] / (double)seededCount;

            fractionFile << t_ << ", "<<seededCount<<", "<<emptyCount<<", ";
            for (int i = 0; i < 16; i++)
                fractionFile << fracs[i] << ", ";

            fractionFile << endl;
        }

        if ((FRAP_t_ > 0) && (FRAP_t_ <= next_grf_t))
        {
            evolve_sync_all_MT_until(FRAP_t_);
            bleach();
            cout << " FRAP at t=" << t_ << "\n";
            FRAP_t_ = -1;
            just_frapped = true;
        }

        if (just_frapped)
        {
            just_frapped = false;
            if (t_ <= next_grf_t - grf_dt_ * 0.01)
                evolve_sync_all_MT_until(next_grf_t);
        }
        else
        {
            evolve_sync_all_MT_until(next_grf_t);
        }

        make_fig(next_grf_t);
        make_histogram(next_grf_t);

        cerr << "t = " << t_ << "\n";
        next_grf_t += grf_dt_;
    }
    if (t_ > next_grf_t - grf_dt_ * 0.5)
    {
        make_fig(next_grf_t);
        make_histogram(next_grf_t);
    };

    cerr << "run done. t=" << t_ << "\n";
    fractionFile.close();
    lengthsFile.close();
    make_kymograph();
}

void Phragmoplast::evolve_sync_all_MT_until(double tmax)
{
    double dt, dt_min;
    int MT_min = -1, event, event_min;
    int trace = 0;

    //int event_no[PHP_EVENT_M_RESEED+1];// MONITOR
    //for(i=0; i<=PHP_EVENT_M_RESEED;i++) { event_no[i] = 0; } // MONITOR

    //std::cerr<<"evolve_sync_all_MT_until("<<tmax<<")\n";
    while (t_ < tmax)
    {
        dt_min = 1e64;
        for (int i = 0; i < N_MT_; ++i)
        {
            event = mt_array_[i].select_next_event(dt);

            if (dt < dt_min)
            {
                dt_min = dt;
                MT_min = i;
                event_min = event;
            }
        }
        if (t_ + dt_min > tmax) // time to stop
        {
            t_ = tmax;
            return;
        }
        //cerr<<"t="<<t_<<" dt="<<dt_min<<" Event="<<event_min<<" i="<<MT_min<<" tot_MT_="<<tot_MT_<<"\n";
        //event_no[event_min]++; // MONITOR
        //cout<<MT_min << ", " << event_min << ", " << tot_MT_0_ << ", " << tot_MT_ << ", " << rel_MT_<<endl;

        mt_array_[MT_min].do_event(event_min);
        t_ += dt_min;

#ifdef MONITOR
        if (trace++ > 10)
        {
            cerr << "t=" << t_ << " last event:" << event_min << "\n";
            trace = 0;
            for (int i = PHP_EVENT_P_GROW; i <= PHP_EVENT_P_RESEED; i++) // MONITOR
            {
                cerr << "event_no[" << i << "]:" << event_no[i] << "\n";
                event_no[i] = 0;
            }
            for (int i = PHP_EVENT_M_GROW; i <= PHP_EVENT_M_RESEED; i++) // MONITOR
            {
                cerr << "event_no[" << i << "]:" << event_no[i] << "\n";
                event_no[i] = 0;
            }
        }
#endif
    }
};

// Bleach each MT one by one
void Phragmoplast::bleach()
{
    for (int i = 0; i < N_MT_; ++i)
    {
        //cerr<<"bleach "<<i;
        mt_array_[i].bleach();
    }
}

// Compute the total MT length and average length
// Write average length in double &average_L
double Phragmoplast::total_MT_length(double &average_L)
{
    double L = 0.0;
    int n = 0;

    for (int i = 0; i < N_MT_; ++i)
    {
        if (mt_array_[i].get_state_p() != PHP_STATE_EMPTY)
        {
            ++n;
            L += mt_array_[i].L();
        }
    }
    average_L = (n > 0) ? L / n : 0;
    return (L);
}

void Phragmoplast::make_fig(double t)
{
    Bitmap bm(Npixel_X_, Npixel_Y_, Thickness_, Width_, t,
              MaxSupperposition_);
    ofstream ofs;
    double empty = 0, mature = 0, growing = 0, shrinking = 0, pause = 0, L, average_L;
    double Grow_out = 0, Grow_in = 0, teadmill = 0;

    //cerr<<"make_fig() Npixel_Y_="<<Npixel_Y_<<"\n";
    for (int i = 0; i < N_MT_; ++i)
    {
        //cerr<<" Draw MT "<<i<<" type="<<mt_array_[i].get_type()
        //<<" state="<<mt_array_[i].get_state_p()<<"\n";
        mt_array_[i].draw(bm);

        if (mt_array_[i].get_state_p() == PHP_STATE_EMPTY)
        {
            ++empty;
        }
        if (mt_array_[i].get_state_p() == PHP_STATE_MATURE)
        {
            ++mature;
        }
        if (mt_array_[i].get_state_p() == PHP_STATE_GROWING)
        {
            ++growing;
        }
        if (mt_array_[i].get_state_p() == PHP_STATE_SHRINKING)
        {
            ++shrinking;
        }
        if (mt_array_[i].get_state_p() == PHP_STATE_PAUSE)
        {
            ++pause;
        }

        if (mt_array_[i].get_type() == PHP_TYPE_GROW_OUTSIDE)
        {
            ++Grow_out;
        }
        if (mt_array_[i].get_type() == PHP_TYPE_GROW_INSIDE)
        {
            ++Grow_in;
        }
        if (mt_array_[i].get_type() == PHP_TYPE_TREADMILL)
        {
            ++teadmill;
        }
    }

    //bm.gaussian_blur(blur_radius_);

    add_kymograph_line(bm);

    //bm.output(fname_.next_name());
    if (!no_grf_)
    {
        bm.save(fname_.next_name());
        save_list(t);
    }
    L = total_MT_length(average_L);

    ofs.open(luminosity_name_, std::ios::app);
    ofs << t << " " << bm.luminosity(Microtuble::FRAPx1__ + Microtuble::dx0__, Microtuble::FRAPy1__ + Microtuble::dx0__, Microtuble::FRAPx2__ - Microtuble::dx0__, Microtuble::FRAPy2__ - Microtuble::dx0__) / N_MT_
        << " " << empty << " " << mature << " " << growing << " " << shrinking << " " << pause << " "
        << Grow_out << " " << Grow_in << " " << teadmill << " "
        << L << " " << average_L << "\n";
    ofs.close();
}

void Phragmoplast::make_histogram(double t)
{
    ofstream ofs;
    double hist[1000];
    int j, N;
    char filename[PHP_PATH_LENGTH + 100];

    //cerr<< "make_histogram("<<t<<")\n";

    sprintf(filename, "%s_t%g.txt", histogram_name_, t);

    for (j = 0; j < bin_N_; ++j)
    {
        hist[j] = 0;
    }
    N = 0;

    for (int i = 0; i < N_MT_; ++i)
    {
        if (mt_array_[i].get_state_p() != PHP_STATE_EMPTY)
        {
            j = (int)floor(mt_array_[i].L() * (bin_N_ - 1) / Microtuble::xmax__);
            if (j >= bin_N_)
                j = bin_N_ - 1;
            hist[j] += 1;
            ++N;
        }
    }

    ofs.open(filename, std::ios::out);
    ofs << "# t=" << t << "\n";
    ofs << "# N=" << N << "\n";
    for (j = 0; j < bin_N_; ++j)
    {
        ofs << j * Microtuble::xmax__ / (bin_N_ - 1.0) << " " << hist[j] / N << "\n";
    }

    ofs.close();
}

void Phragmoplast::save_list(double t)
{
    ofstream ofs(fname_.current_name_other_suffix(".txt"), std::ios::app);

    ofs << "#index type state x1 y1 x2 y2 bleached_x1 bleached_x2 \n";
    ofs << "#t=" << t << "\n";
    for (int i = 0; i < N_MT_; ++i)
    {
        ofs << i << " ";
        mt_array_[i].txt_save(ofs);
    }
}

void Phragmoplast::add_kymograph_line(Bitmap &bm)
{
    double v[Npixel_X_];

    bm.average(v, kymograph_amp_);

    kymograph_bm_->add_line(v, index_kymograph_++);
}

void Phragmoplast::make_kymograph()
{
    kymograph_bm_->output(kymograph_name_);
}

double Phragmoplast::R_gs(double L_MT, double x)
{
    double r_gs = r_gs_ / (rel_MT_ + 1e-6);

    if (L_MT > macet_MT_len_)
    {
        double macet_f = macet_mz_ + (macet_mz_ - macet_dz_) *
                                         (x - Width_ + macet_MT_len_) / macet_MT_len_;
        r_gs *= 1 + macet_f * (L_MT - macet_MT_len_) / macet_MT_len_; //possibly the Width_
    }
    // if(L_MT < .2)
    // {
    //     r_gs *= .1;
    // }
    return (r_gs);
}

double Phragmoplast::R_gp(double L_MT, double x)
{
    double r_gp = r_gp_ / (rel_MT_ + 1e-6);

    if (L_MT > macet_MT_len_)
    {
        double macet_f = macet_mz_ + (macet_mz_ - macet_dz_) *
                                         (x - Width_ + macet_MT_len_) / macet_MT_len_;
        r_gp *= 1 + macet_f * (L_MT - macet_MT_len_) / macet_MT_len_;
    }
    // if(L_MT < .2)
    // {
    //     r_gp *= .1;
    // }
    return (r_gp);
}

double Phragmoplast::R_me_gs(double L_MT, double x)
{
    double r_me_gs = r_me_gs_ / (rel_MT_ + 1e-6);

    if (L_MT > macet_MT_len_)
    {
        double macet_f = macet_mz_ + (macet_mz_ - macet_dz_) *
                                         (x - Width_ + macet_MT_len_) / macet_MT_len_;
        r_me_gs *= 1 + macet_f * (L_MT - macet_MT_len_) / macet_MT_len_;
    }
    // if(L_MT < .2)
    // {
    //     r_me_gs *= .1;
    // }
    return (r_me_gs);
}

double Phragmoplast::R_me_gp(double L_MT, double x)
{
    double r_me_gp = r_me_gp_ / (rel_MT_ + 1e-6);

    if (L_MT > macet_MT_len_)
    {
        double macet_f = macet_mz_ + (macet_mz_ - macet_dz_) *
                                         (x - Width_ + macet_MT_len_) / macet_MT_len_;
        r_me_gp *= 1 + macet_f * (L_MT - macet_MT_len_) / macet_MT_len_;
    }
    // if(L_MT < .2)
    // {
    //     r_me_gp *= .1;
    // }
    return (r_me_gp);
}
