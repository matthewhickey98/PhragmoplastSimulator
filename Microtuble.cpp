#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "Microtuble.h"
#include "Phragmoplast.h"
#include <Magick++.h>
using namespace std;
//using namespace Magick;

// need to do this for them to be seen.
double Microtuble::xmax__, Microtuble::ymax__, Microtuble::dx0__, Microtuble::min_L__;
double Microtuble::FRAPx1__, Microtuble::FRAPy1__, Microtuble::FRAPx2__, Microtuble::FRAPy2__;
bool Microtuble::NoMature_mode__;

//************************************************************/
// class mt
void Microtuble::init(double x1, double y1, double L, double theta, int type,
                      Phragmoplast *php)
{
    php_ = php;

    type_ = type;
    state_p_ = PHP_STATE_GROWING;
    state_m_ = PHP_STATE_EMPTY;
    towards_DZ_ = type == PHP_TYPE_GROW_OUTSIDE;
    if (type_ == PHP_TYPE_TREADMILL)
    {
        state_m_ = PHP_STATE_GROWING;
    }

    if (towards_DZ_)
    {
        x1_ = x1;
        y1_ = y1;

        x2_ = x1_ + L * cos(theta * M_PI / 180.0);
        y2_ = y1_ + L * sin(theta * M_PI / 180.0);
    }
    else
    {
        x2_ = x1;
        y2_ = y1;

        x1_ = x2_ - L * cos(theta * M_PI / 180.0);
        y1_ = y2_ - L * sin(theta * M_PI / 180.0);
    }

    dx_ = dx0__ * cos(theta * M_PI / 180.0);
    dy_ = dx0__ * sin(theta * M_PI / 180.0);

    //  cerr<<"theta="<<theta<<" L = "<<L<<" x1="<<x1_<<" y1="<<y1_<<" x2_="<<x2_<<" y2_="<<y2_<<" tan(theta)="<<(y2_-y1_)/(x2_-x1)<<"\n";
    check_mt_pos();
}

double Microtuble::rnd(double x)
{
    return (php_->rnd(x));
}

// Check MT pos if necessary and trim it if needed.
void Microtuble::check_mt_pos()
{
   if(x1_ < 0) 
  { if(x2_ > Microtuble::min_L__)
    { y1_ = line_y(0);
      if (! towards_DZ_) // NEW
      { state_p_ = PHP_STATE_IN_CP;
      }
    } 
    else 
    { x2_ = 0; 
      state_p_ = PHP_STATE_EMPTY; 
    } 
    x1_ = 0; 
  } 
  if(x1_ > Microtuble::xmax__) 
  { x1_ = x2_ = Microtuble::xmax__; 
    state_p_ = PHP_STATE_EMPTY;  
  }

  if(x2_ < 0) 
  { x1_ = x2_ = 0; 
    state_p_ = PHP_STATE_EMPTY; 
  } 
  if(x2_ > Microtuble::xmax__) 
  { if(x1_ < (Microtuble::xmax__-Microtuble::min_L__))
    { y2_ = line_y(Microtuble::xmax__);
      if (towards_DZ_) // NEW
      { state_p_ = PHP_STATE_IN_DZ;
      }

    }
    else 
    { x1_ = Microtuble::xmax__; 
      state_p_ = PHP_STATE_EMPTY; 
    } 
    x2_ = Microtuble::xmax__; 
  }

  if(y1_ < 0) 
  { if(y2_ > Microtuble::min_L__) 
    { x1_ = line_x(0); 
      if(!NoMature_mode__) { state_p_ = PHP_STATE_MATURE; } 
    } 
    else { y2_ = 0;  } 
    y1_ = 0; 
  } 
  if(y1_ > Microtuble::ymax__) 
  { if(y2_ < (Microtuble::ymax__-Microtuble::min_L__)) 
    { x1_ = line_x(Microtuble::ymax__); 
      if(!NoMature_mode__) { state_p_ = PHP_STATE_MATURE; }    
    } 
    else { y2_ = Microtuble::ymax__;  } 
    y1_ = Microtuble::ymax__; 
  }

  if(y2_ < 0) 
  { if(y1_ > Microtuble::min_L__) 
    { x2_ = line_x(0); 
     if(!NoMature_mode__) { state_p_ = PHP_STATE_MATURE; } 
    } 
    else 
    { y1_ = 0; 
    } 
    y2_ = 0;   
  } 
  if(y2_ > Microtuble::ymax__) 
  { if(y1_ < (Microtuble::ymax__-Microtuble::min_L__)) 
    { x2_ = line_x(Microtuble::ymax__); 
      if(!NoMature_mode__) { state_p_ = PHP_STATE_MATURE; }
    } 
    else 
    { y1_ = Microtuble::ymax__;  
    } 
    y2_ = Microtuble::ymax__;
  }
 
  if(state_p_ != PHP_STATE_EMPTY) 
  { if(  (x1_ >= (x2_-Microtuble::min_L__))
	 ||((x2_ < x1_) || (L() < Microtuble::min_L__))
	 ||(((x1_ > (Microtuble::xmax__-Microtuble::min_L__))
	     && (x2_ > (Microtuble::xmax__-Microtuble::min_L__)))
            ||((x1_ < Microtuble::min_L__) && (x2_ < Microtuble::min_L__))))
    { if(towards_DZ_)  { x2_ = x1_; }
       else             { x1_ = x2_; }
       state_p_ = PHP_STATE_EMPTY;
       //cerr <<" set PHP_STATE_EMPTY 1"<<"state="<<state_p_<<"\n";
    }
  }

  if(state_p_ == PHP_STATE_EMPTY)
  { bleached_x1_ = bleached_x2_ = -1;
  }

}

void Microtuble::draw(Bitmap &bm)
{
    //cerr<<"draw x1="<<x1_<<" y1="<<y1_<<" x2="<<x2_<<" y2="<<y2_
    //    <<" bleached_x1_="<<bleached_x1_<<" bleached_x2_="<<bleached_x2_
    //    <<" type="<<type_<<" state="<<state_p_<<"\n";

    if (state_p_ == PHP_STATE_EMPTY)
    {
        return;
    }

    if (bleached_x1_ >= 0) // MT is bleached
    {
        if (x1_ < bleached_x1_)
        {
            bm.line(x1_, y1_, bleached_x1_, line_y(bleached_x1_));
        }
        if (x2_ > bleached_x2_)
        {
            bm.line(bleached_x2_, line_y(bleached_x2_), x2_, y2_);
        }
    }
    else // not bleached at all
    {
        bm.line(x1_, y1_, x2_, y2_);
    }
};

void Microtuble::update_bleach()
{

    //cerr<<"update_bleach("<<dx<<")  type="<<type_<<"\n";
    if ((bleached_x1_ >= 0) && (bleached_x2_ >= 0))
    {
        if (bleached_x2_ > x2_)
        {
            bleached_x2_ = x2_;
        }
        if (bleached_x1_ > x2_)
        {
            bleached_x1_ = bleached_x2_ = -1;
        }
        if (bleached_x1_ < x1_)
        {
            bleached_x1_ = x1_;
        }
        if (bleached_x2_ < x1_)
        {
            bleached_x1_ = bleached_x2_ = -1;
        }
    }
}

// Adjust the size of theMT
// plus true : increase size by dx,dy at + end
//      false: increase size by dx,dy at - end
// grow true : increase size by dx,dy
//      false: decrease size by dx,dy
// Test MT bleach and position.
void Microtuble::size_adjust(bool plus, bool grow)
{

    //std::cerr<<"towards_DZ_="<<towards_DZ_<<" plus="<<plus<<" grow="<<grow<<"\n";

    if (grow)
    {
        if (towards_DZ_ == plus)
        {
            x2_ += dx_;
            y2_ += dy_;
        } // Growing from right;
        else
        {
            x1_ -= dx_;
            y1_ -= dy_;
            //std::cerr<<"Growing from left  x1_="<<x1_<<" x2_="<<x2_<<" dx_="<<dx_<<" dy_="<<dy_<<"\n";
        } // Growing from left;
        php_->tot_MT_ -= dx0__;
    }
    else
    {
        if (towards_DZ_ == plus)
        {
            x2_ -= dx_;
            y2_ -= dy_;
        } // Shrinking from right;
        else
        {
            x1_ += dx_;
            y1_ += dy_;
            //std::cerr<<"Shrinking from left\n";
        } // Shrinking from left;
        php_->tot_MT_ += dx0__;
    }

    php_->rel_MT_ = php_->tot_MT_ / php_->tot_MT_0_;

    // test bleach reagion
    if ((bleached_x1_ >= 0) && (bleached_x2_ >= 0))
    {
        if (bleached_x2_ > x2_)
        {
            bleached_x2_ = x2_;
        }
        if (bleached_x1_ > x2_)
        {
            bleached_x1_ = bleached_x2_ = -1;
        }
        if (bleached_x1_ < x1_)
        {
            bleached_x1_ = x1_;
        }
        if (bleached_x2_ < x1_)
        {
            bleached_x1_ = bleached_x2_ = -1;
        }
    }

    check_mt_pos();
}

void Microtuble::do_event(int event)
{
    if (event & PHP_EVENT_MASK_MINUS) // minus end event
    {
        switch (state_m_)
        {
        case PHP_STATE_GROWING:
            switch (event)
            {
            case PHP_EVENT_M_GROW:
                if (dx0__ < php_->tot_MT_)
                {
                    size_adjust(false, true);
                }
                return;

            case PHP_EVENT_M_CATASTROPHE:
                state_m_ = PHP_STATE_SHRINKING;
                size_adjust(false, false);
                return;

            case PHP_EVENT_M_PAUSE:
                state_m_ = PHP_STATE_PAUSE;
                return;
            }
            break;

        case PHP_STATE_SHRINKING:
            switch (event)
            {
            case PHP_EVENT_M_SHORTEN:
                size_adjust(false, false);
                return;

            case PHP_EVENT_M_RECOVERY:
                if (dx0__ <= php_->tot_MT_)
                {
                    state_m_ = PHP_STATE_GROWING; // NEW: move this inside test
                    size_adjust(false, true);
                }
                return;

            case PHP_EVENT_M_PAUSE:
                state_m_ = PHP_STATE_PAUSE;
                return;
            }
            break;

        case PHP_STATE_PAUSE:
            switch (event)
            {
            case PHP_EVENT_M_RECOVERY:
                if (dx0__ <= php_->tot_MT_)
                {
                    state_m_ = PHP_STATE_GROWING; // NEW: move this inside test
                    size_adjust(false, true);
                }
                return;

            case PHP_EVENT_M_CATASTROPHE:
                state_m_ = PHP_STATE_SHRINKING;
                size_adjust(false, false);
                return;
            }
            break;

        case PHP_STATE_MATURE:
            if (event == PHP_EVENT_M_CATASTROPHE)
            {
                state_m_ = PHP_STATE_SHRINKING;
                size_adjust(false, false);
            }
            return;
        }
    }
    else // event at + end
    {
        switch (state_p_)
        {
        case PHP_STATE_EMPTY:
            if ((event == PHP_EVENT_P_RESEED) && (dx0__ < php_->tot_MT_))
            {
                theta_ = (php_->rnd(1.0) - 0.5) * 2.0 * php_->get_theta_max();
                x1_ = x2_ = php_->rnd_.next() * php_->get_Thickness();
                y1_ = y2_ = php_->rnd_.next() * php_->get_Width();
                dx_ = dx0__ * cos(theta_ * M_PI / 180.0);
                dy_ = dx0__ * sin(theta_ * M_PI / 180.0);
                if (towards_DZ_)
                {
                    x2_ += dx_;
                    y2_ += dy_;
                }
                else
                {
                    x1_ -= dx_;
                    y1_ -= dy_;
                }

                size_adjust(true, true);
                state_p_ = PHP_STATE_GROWING;
            }
            return;

        case PHP_STATE_GROWING:
            switch (event)
            {
            case PHP_EVENT_P_GROW:
                if (dx0__ < php_->tot_MT_)
                {
                    size_adjust(true, true);
                }
                return;

            case PHP_EVENT_P_CATASTROPHE:
                state_p_ = PHP_STATE_SHRINKING;
                size_adjust(true, false);
                return;

            case PHP_EVENT_P_PAUSE:
                state_p_ = PHP_STATE_PAUSE;
                return;
            }
            break;

        case PHP_STATE_SHRINKING:
            switch (event)
            {
            case PHP_EVENT_P_SHORTEN:
                size_adjust(true, false);
                return;

            case PHP_EVENT_P_RECOVERY:
                if (dx0__ <= php_->tot_MT_)
                {
                    state_p_ = PHP_STATE_GROWING; // NEW: move this inside test
                    size_adjust(true, true);
                }
                return;

            case PHP_EVENT_P_PAUSE:
                state_p_ = PHP_STATE_PAUSE;
                return;
            }
            break;

        case PHP_STATE_PAUSE:
            switch (event)
            {
            case PHP_EVENT_P_RECOVERY:
                if (dx0__ <= php_->tot_MT_)
                {
                    state_p_ = PHP_STATE_GROWING; // NEW: move this inside test
                    size_adjust(true, true);
                }
                return;

            case PHP_EVENT_P_CATASTROPHE:
                state_p_ = PHP_STATE_SHRINKING;
                size_adjust(true, false);
                return;
            }
            break;

        case PHP_STATE_MATURE:
            if (event == PHP_EVENT_P_CATASTROPHE)
            {
                state_p_ = PHP_STATE_SHRINKING;
                size_adjust(true, false);
            }
            return;

        case PHP_STATE_IN_CP:
        case PHP_STATE_IN_DZ:
            if (event == PHP_EVENT_P_CATASTROPHE)
            {
                state_p_ = PHP_STATE_SHRINKING;
                size_adjust(true, false);
            }
            return;
        }
    }
}

#define EVENT(RATE, TYPE)           \
    {                               \
        dtt = php_->random_t(RATE); \
        if (dtt < dt_min)           \
        {                           \
            event_min = TYPE;       \
            dt_min = dtt;           \
        }                           \
    }

// select which event to do next (Monte Carlo) for this MT
// Write time into dt and return event to implement
int Microtuble::select_next_event(double &dt)
{
    int event_min = 0;
    double dtt, dt_min = 1e64, L_MT, xp, xm;

    //  cerr<<"state ="<<state_p_<<":\n";

    // Polymerisation of + end
    dt_min = 1e32;
    event_min = PHP_EVENT_NONE;

    L_MT = L(); // length of this MT
    if (towards_DZ_)
    {
        xp = x2_;
        xm = x1_;
    }
    else
    {
        xm = x2_;
        xp = x1_;
    }

    // Polymerisation of - end
    if (type_ == PHP_TYPE_TREADMILL)
    {
        switch (state_m_)
        {
        case PHP_STATE_GROWING:
            EVENT(php_->R_me_polym() / dx0__, PHP_EVENT_M_GROW);
            EVENT(php_->R_me_gs(L_MT, xm), PHP_EVENT_M_CATASTROPHE);
            EVENT(php_->R_me_gp(L_MT, xm), PHP_EVENT_M_PAUSE);
            break;
        case PHP_STATE_SHRINKING:
            EVENT(php_->get_r_me_depolym() / dx0__, PHP_EVENT_M_SHORTEN);
            EVENT(php_->R_me_sg(), PHP_EVENT_M_RECOVERY);
            EVENT(php_->R_me_sp(), PHP_EVENT_M_PAUSE);
            break;

        case PHP_STATE_PAUSE:
            EVENT(php_->R_me_pg(), PHP_EVENT_M_RECOVERY);
            EVENT(php_->R_me_ps(), PHP_EVENT_M_CATASTROPHE);
            break;

        case PHP_STATE_MATURE:
            EVENT(php_->R_me_gs(L_MT, xm), PHP_EVENT_M_CATASTROPHE);
            break;
        default:
            cerr << "invalid state : " << state_p_ << "\n";
            exit(0);
        }
        dt = dt_min;
    }

    // + end
    switch (state_p_)
    {
    case PHP_STATE_EMPTY: // to improve
        EVENT(php_->R_r_reseed(), PHP_EVENT_P_RESEED);
        break;

    case PHP_STATE_GROWING:
        EVENT(php_->R_polym() / dx0__, PHP_EVENT_P_GROW);
        EVENT(php_->R_gs(L_MT, xp), PHP_EVENT_P_CATASTROPHE);
        EVENT(php_->R_gp(L_MT, xp), PHP_EVENT_P_PAUSE);
        break;
    case PHP_STATE_SHRINKING:
        EVENT(php_->get_r_depolym()/ dx0__, PHP_EVENT_P_SHORTEN); //was get_r_me_depolym()
        EVENT(php_->R_sg(), PHP_EVENT_P_RECOVERY);
        EVENT(php_->R_sp(), PHP_EVENT_P_PAUSE);
        break;

    case PHP_STATE_PAUSE:
        EVENT(php_->R_pg(), PHP_EVENT_P_RECOVERY);
        EVENT(php_->R_ps(), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_MATURE:
        EVENT(php_->R_gs(L_MT, xp), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_IN_CP:
        EVENT(php_->R_ps_CP(), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_IN_DZ:
        EVENT(php_->R_ps_DZ(), PHP_EVENT_P_CATASTROPHE);
        break;

    default:
        cerr << "invalid state : " << state_p_ << "\n";
        exit(0);
    }
    dt = dt_min;

    //  cerr <<"event: "<<event_min<< " dt="<<dt_min<<"\n";
    return (event_min);
}

// Evolve this MT until t_max
void Microtuble::run_until(double t_max)
{
    double dt;
    int event;

    while (t_ < t_max)
    {
        event = select_next_event(dt);
        //cerr<<" t="<<t_<<" event="<<event<<" dt="<<dt<<"\n";
        //    if(t_+dt >= t_max) { return; }
        do_event(event);
        t_ += dt;
    }
}

// y = ((y2-y1)*x+y1*x2-y2*x1)/(x2-x1)
// x = ((x2-x1)*y+x1*y2-x2*y1)/(y2-y1)
void Microtuble::bleach()
{
    double x, y;
    int i = 0;

    //cerr<<" x1="<<x1_<<" y1="<<y1_<<" x2="<<x2_<<" y2="<<y2_<<" type="<<type_<<"\n";

    if (type_ == PHP_TYPE_GROW_INSIDE)
    {
        // swap x1_ and x2_ temporaily
        x = x1_;
        x1_ = x2_;
        x2_ = x;
        y = y1_;
        y1_ = y2_;
        y2_ = y;
    }

    //cerr<<"Bleach x1=("<<x1_<<","<<y1_<<")  x2=("<<x2_<<","<<y2_<<")\n";

    if (is_inside_FRAP(x1_, y1_))
    {
        bleached_x1_ = x1_;
        ++i;
        //cerr<<"x1,y1 Inside\n";
    }
    else
    {
        //cerr<<"x1,y1 Not Inside\n";
    }

    if (is_inside_FRAP(x2_, y2_))
    {
        if (i == 0)
        {
            bleached_x1_ = x2_;
        }
        else
        {
            bleached_x2_ = x2_;
        }
        ++i;
        //cerr<<"x2,y2 Inside\n";
    }
    else
    {
        //cerr<<"x2,y2 Not Inside\n";
    }

    if (i < 2)
    {
        if ((y = cross_left_FRAP()) >= 0)
        {
            if (i == 0)
            {
                bleached_x1_ = FRAPx1__;
            }
            else
            {
                bleached_x2_ = FRAPx1__;
            }
            ++i;
            //cerr<<"cross_left_FRAP y="<<y<<"\n";
        }

        if ((y = cross_right_FRAP()) >= 0)
        {
            if (i == 0)
            {
                bleached_x1_ = FRAPx2__;
            }
            else
            {
                bleached_x2_ = FRAPx2__;
            }
            ++i;
            //cerr<<"cross_right_FRAP y="<<y<<"\n";
        }

        if ((x = cross_bottom_FRAP()) >= 0)
        {
            if (i == 0)
            {
                bleached_x1_ = x;
            }
            else
            {
                bleached_x2_ = x;
            }
            ++i;
            //cerr<<"cross_bottom_FRAP y="<<x<<"\n";
        }

        if ((x = cross_top_FRAP()) >= 0)
        {
            if (i == 0)
            {
                bleached_x1_ = x;
            }
            else
            {
                bleached_x2_ = x;
            }
            ++i;
            //cerr<<"cross_top_FRAP y="<<x<<"\n";
        }
    }

    if (i >= 2)
    {
        if (bleached_x1_ > bleached_x2_)
        {
            x = bleached_x1_;
            bleached_x1_ = bleached_x2_;
            bleached_x2_ = x;
        }
    }
    else
    {
        bleached_x1_ = bleached_x2_ = -1;
        //cerr <<"Unbleached x1="<<x1_<<", y1="<<y1_<<" x2="<<x2_<<", y2="<<y2_<<"\n";
    }

    //cerr<<"  bleach x1="<<bleached_x1_<<" x2="<<bleached_x2_<<"\n";

    if (type_ == PHP_TYPE_GROW_INSIDE)
    {
        // swap x1_ and x2_  back
        x = x1_;
        x1_ = x2_;
        x2_ = x;
        y = y1_;
        y1_ = y2_;
        y2_ = y;
    }
}

// Save MT as 1 line in txt file
// type state x1 y1 x2 y2 bx1 bx2
void Microtuble::txt_save(ofstream &ofs)
{
    ofs << type_ << " " << state_p_ << " "
        << x1_ << " " << y1_ << " " << x2_ << " " << y2_ << " "
        << bleached_x1_ << " " << bleached_x2_ << "\n";
}
