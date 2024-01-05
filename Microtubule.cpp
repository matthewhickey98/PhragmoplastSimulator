#include "Microtubule.hpp"
// #include <Magick++.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <string>

using namespace std;

/**
 * @brief Return the time of an event given the rate and type of event
 *
 */
#define EVENT(RATE, TYPE)                     \
    {                                         \
        dtt = phragmoplast->RandomTime(RATE); \
        if (dtt < dt_min)                     \
        {                                     \
            event_min = TYPE;                 \
            dt_min = dtt;                     \
        }                                     \
    }

// using namespace Magick;

//************************************************************/
// class mt

/**
 * @brief Initialize the microtubule
 *
 * @param x1        x-coordinate of the seed
 * @param y1        y-coordinate of the seed
 * @param L         length of the microtubule
 * @param theta     angle of the microtubule
 * @param type      type of the microtubule
 * @param php       pointer to the phragmoplast
 */
void Microtubule::Init(double x1, double y1, double L, double theta, int type, Phragmoplast *php)
{
    phragmoplast = php;
    xInit = x1;
    yInit = y1;
    typeMT = type;

    statePlusEnd = PHP_STATE_GROWING;
    stateMinusEnd = PHP_STATE_EMPTY;

    towardsDistalZone = type == PHP_TYPE_GROW_OUTSIDE;

    if (typeMT == PHP_TYPE_TREADMILL)
    {
        stateMinusEnd = PHP_STATE_GROWING;
    }

    if (towardsDistalZone)
    {
        this->x1 = x1;
        this->y1 = y1;

        this->x2 = x1 + L * cos(theta * M_PI / 180.0);
        this->y2 = y1 + L * sin(theta * M_PI / 180.0);
    }
    else
    {
        this->x2 = x1;
        this->y2 = y1;

        this->x1 = x2 - L * cos(theta * M_PI / 180.0);
        this->y1 = y2 - L * sin(theta * M_PI / 180.0);
    }

    deltaX = phragmoplast->DeltaX * cos(theta * M_PI / 180.0);
    deltaY = phragmoplast->DeltaX * sin(theta * M_PI / 180.0);

    CheckMTPosition();
}

/**
 * @brief Check the microtubule position and trim it if necessary
 *
 */
void Microtubule::CheckMTPosition()
{
    // bleachedx1 = bleachedx2 = -1;
    if (statePlusEnd == PHP_STATE_EMPTY)
    {
        x1 = x2 = xInit;
        y1 = y2 = yInit;
        bleachedx1 = bleachedx2 = -1;
        return;
    }
    if ((x1 <= phragmoplast->CellPlateXMax && typeMT == PHP_TYPE_GROW_INSIDE)) //|| x2 <= CP_X_MAX)
    {
        statePlusEnd = PHP_STATE_IN_CP;
    }
    else if ((x2 >= phragmoplast->DistalZoneXMin && typeMT == PHP_TYPE_GROW_OUTSIDE))
    {
        statePlusEnd = PHP_STATE_IN_DZ;
    }

    if (x1 < 0)
    {
        if (x2 > phragmoplast->DeltaX * 0.05)
        {
            y1 = LineY(0);
        }
        else
        {
            x2 = 0;
            statePlusEnd = PHP_STATE_EMPTY;
        }
        x1 = 0;
    }
    if (x1 > phragmoplast->ThicknessMT)
    {
        x1 = phragmoplast->ThicknessMT;
        statePlusEnd = PHP_STATE_IN_DZ;
    }

    if (x2 < 0)
    {
        x2 = 0;
        // state_p_ = PHP_STATE_IN_CP;
    }
    if (x2 > phragmoplast->ThicknessMT)
    {
        if (x1 < (phragmoplast->ThicknessMT - phragmoplast->DeltaX * 0.05))
        {
            y2 = LineY(phragmoplast->ThicknessMT);
            statePlusEnd = PHP_STATE_IN_DZ;
        }
        x2 = phragmoplast->ThicknessMT;
    }

    if (y1 < 0)
    {
        if (y2 > phragmoplast->DeltaX * 0.05)
        {
            x1 = LineX(0);
            if (!phragmoplast->NoMatureMode)
            {
                statePlusEnd = PHP_STATE_MATURE;
            }
        }

        y1 = 0;
    }
    if (y1 > phragmoplast->WidthMT)
    {
        if (y2 < (phragmoplast->WidthMT - phragmoplast->DeltaX * 0.05))
        {
            x1 = LineX(phragmoplast->WidthMT);
            if (!phragmoplast->NoMatureMode)
            {
                statePlusEnd = PHP_STATE_MATURE;
            }
        }
        else
        {
            y2 = phragmoplast->WidthMT;
            statePlusEnd = PHP_STATE_EMPTY;
        } // 2MAR21
        y1 = phragmoplast->WidthMT;
    }

    if (y2 < 0)
    {
        if (y1 > phragmoplast->DeltaX * 0.05)
        {
            x2 = LineX(0);
            if (!phragmoplast->NoMatureMode)
            {
                statePlusEnd = PHP_STATE_MATURE;
            }
        }
        else
        {
            y1 = 0;
            statePlusEnd = PHP_STATE_EMPTY; // 2MAR21
        }
        y2 = 0;
    }
    if (y2 > phragmoplast->WidthMT)
    {
        if (y1 < (phragmoplast->WidthMT - phragmoplast->DeltaX * 0.05))
        {
            x2 = LineX(phragmoplast->WidthMT);
            if (!phragmoplast->NoMatureMode)
            {
                statePlusEnd = PHP_STATE_MATURE;
            }
        }
        else
        {
            y1 = phragmoplast->WidthMT;
            statePlusEnd = PHP_STATE_EMPTY; // 2MAR21
        }
        y2 = phragmoplast->WidthMT;
    }

    if (statePlusEnd != PHP_STATE_EMPTY)
    {
        if ((x1 >= (x2 - phragmoplast->DeltaX * 0.05)) || ((x2 < x1) || (Length() < phragmoplast->DeltaX * 0.05)) ||
            (((x1 > (phragmoplast->ThicknessMT - phragmoplast->DeltaX * 0.05)) &&
              (x2 > (phragmoplast->ThicknessMT - phragmoplast->DeltaX * 0.05))) ||
             ((x1 < phragmoplast->DeltaX * 0.05) && (x2 < phragmoplast->DeltaX * 0.05))))
        {
            if (towardsDistalZone)
            {
                x2 = x1;
            }
            else
            {
                x1 = x2;
            }
            statePlusEnd = PHP_STATE_EMPTY;
            // cout <<" set PHP_STATE_EMPTY 1"<<"state="<<state_p_<<"\n";
        }
    }
    if ((bleachedx1 >= 0) && (bleachedx2 >= 0))
    {
        if (bleachedx2 > x2)
        {
            bleachedx2 = x2;
        }
        if (bleachedx1 > x2)
        {
            bleachedx1 = bleachedx2 = -1;
        }
        if (bleachedx1 < x1)
        {
            bleachedx1 = x1;
        }
        if (bleachedx2 < x1)
        {
            bleachedx1 = bleachedx2 = -1;
        }
    }
    if (statePlusEnd == PHP_STATE_EMPTY)
    {
        bleachedx1 = bleachedx2 = -1;
    }
}

/**
 * @brief Draw a microtubule on the bitmap
 *
 * @param bm    Bitmap to draw on
 */
void Microtubule::Draw(Bitmap &bm)
{
    if (statePlusEnd == PHP_STATE_EMPTY)
    {
        return;
    }

    if (bleachedx1 >= 0) // MT is bleached
    {
        if (x1 < bleachedx1)
        {
            bm.Line(x1, y1, bleachedx1, LineY(bleachedx1));
        }
        if (x2 > bleachedx2)
        {
            bm.Line(bleachedx2, LineY(bleachedx2), x2, y2);
        }
    }
    else // not bleached at all
    {
        bm.Line(x1, y1, x2, y2);
    }
};

/**
 * @brief Update the coordinates of the microtubule
 *
 */
void Microtubule::UpdateBleach()
{
    if ((bleachedx1 >= 0) && (bleachedx2 >= 0))
    {
        if (bleachedx2 > x2)
        {
            bleachedx2 = x2;
        }
        if (bleachedx1 > x2)
        {
            bleachedx1 = bleachedx2 = -1;
        }
        if (bleachedx1 < x1)
        {
            bleachedx1 = x1;
        }
        if (bleachedx2 < x1)
        {
            bleachedx1 = bleachedx2 = -1;
        }
    }
}

// Adjust the size of theMT
// plus true : increase size by dx,dy at + end
//      false: increase size by dx,dy at - end
// grow true : increase size by dx,dy
//      false: decrease size by dx,dy
// Test MT Bleach and position.

/**
 * @brief Adjust the size of the microtubule
 *
 * @param plus  true : increase size by dx,dy at + end, false: increase size by dx,dy at - end
 * @param grow  true : increase size by dx,dy, false: decrease size by dx,dy
 */
void Microtubule::SizeAdjust(bool plus, bool grow)
{
    double Lold = Length();

    if (grow)
    {
        if (towardsDistalZone == plus)
        {
            x2 += deltaX;
            y2 += deltaY;
        }
        else
        {
            x1 -= deltaX;
            y1 -= deltaY;
        }
    }
    else
    {
        if (towardsDistalZone == plus)
        {
            x2 -= deltaX;
            y2 -= deltaY;
        }
        else
        {
            x1 += deltaX;
            y1 += deltaY;
        }
    }

    // test Bleach region
    if ((bleachedx1 >= 0) && (bleachedx2 >= 0))
    {
        if (bleachedx2 > x2)
        {
            bleachedx2 = x2;
        }
        if (bleachedx1 > x2)
        {
            bleachedx1 = bleachedx2 = -1;
        }
        if (bleachedx1 < x1)
        {
            bleachedx1 = x1;
        }
        if (bleachedx2 < x1)
        {
            bleachedx1 = bleachedx2 = -1;
        }
    }

    CheckMTPosition();
    phragmoplast->IncrementTotalLength(Lold - Length());
    phragmoplast->RelativeMT();
}

/**
 * @brief Do an event on the microtubule
 *
 * @param event     Event to do
 * @return int      0: error, 1: success
 */
int Microtubule::DoEvent(int event)
{
    if (event & PHP_EVENT_MASK_MINUS) // minus end event
    {
        switch (stateMinusEnd)
        {
        case PHP_STATE_GROWING:
            switch (event)
            {
            case PHP_EVENT_M_GROW:
                if (phragmoplast->DeltaX < phragmoplast->GetTotalLength())
                {
                    SizeAdjust(false, true);
                }
                return 1;

            case PHP_EVENT_M_CATASTROPHE:
                stateMinusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(false, false);
                return 1;

            case PHP_EVENT_M_PAUSE:
                stateMinusEnd = PHP_STATE_PAUSE;
                return 1;
            }
            break;

        case PHP_STATE_SHRINKING:
            switch (event)
            {
            case PHP_EVENT_M_SHORTEN:
                SizeAdjust(false, false);
                return 1;

            case PHP_EVENT_M_RECOVERY:
                if (phragmoplast->DeltaX <= phragmoplast->GetTotalLength())
                {
                    stateMinusEnd = PHP_STATE_GROWING;
                    SizeAdjust(false, true);
                }
                return 1;

            case PHP_EVENT_M_PAUSE:
                stateMinusEnd = PHP_STATE_PAUSE;
                return 1;
            }
            break;

        case PHP_STATE_PAUSE:
            switch (event)
            {
            case PHP_EVENT_M_RECOVERY:
                if (phragmoplast->DeltaX <= phragmoplast->GetTotalLength())
                {
                    stateMinusEnd = PHP_STATE_GROWING;
                    SizeAdjust(false, true);
                }
                return 1;

            case PHP_EVENT_M_CATASTROPHE:
                stateMinusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(false, false);
                return 1;
            }
            break;

        case PHP_STATE_MATURE:
            if (event == PHP_EVENT_M_CATASTROPHE)
            {
                stateMinusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(false, false);
            }
            return 1;
        }
    }
    else // event at + end
    {
        switch (statePlusEnd)
        {
        case PHP_STATE_EMPTY:
            if ((event == PHP_EVENT_P_RESEED) && (phragmoplast->DeltaX < phragmoplast->GetTotalLength()))
            {
                theta = (phragmoplast->Random(1.0) - 0.5) * 2.0 * phragmoplast->GetThetaMax();

                // CHANGED to always have the same seed location and to reseed at random length [0, L0]
                // to be consistant with how MT are seeded/reseeded
                // this was necesary to reproduce length distribution
                x1 = x2 = xInit;
                y1 = y2 = yInit;
                double L = phragmoplast->Next() * phragmoplast->InitialLength;
                if (typeMT == PHP_TYPE_GROW_OUTSIDE)
                {
                    x1 = xInit;
                    y1 = yInit;

                    x2 = x1 + L * cos(theta * M_PI / 180.0);
                    y2 = y1 + L * sin(theta * M_PI / 180.0);
                }
                else
                {
                    x2 = xInit;
                    y2 = yInit;

                    x1 = x2 - L * cos(theta * M_PI / 180.0);
                    y1 = y2 - L * sin(theta * M_PI / 180.0);
                }
                phragmoplast->DecrementTotalLength(L);
                bleachedx1 = bleachedx2 = -1;
                deltaX = phragmoplast->DeltaX * cos(theta * M_PI / 180.0);
                deltaY = phragmoplast->DeltaX * sin(theta * M_PI / 180.0);
                statePlusEnd = PHP_STATE_GROWING;
                SizeAdjust(towardsDistalZone, true); // 2MAR21
            }
            return 1;

        case PHP_STATE_GROWING:
            switch (event)
            {
            case PHP_EVENT_P_GROW:
                if (phragmoplast->DeltaX < phragmoplast->GetTotalLength())
                {
                    SizeAdjust(true, true);
                }
                return 1;

            case PHP_EVENT_P_CATASTROPHE:
                statePlusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(true, false);
                return 1;

            case PHP_EVENT_P_PAUSE:
                statePlusEnd = PHP_STATE_PAUSE;
                return 1;
            }
            break;

        case PHP_STATE_SHRINKING:
            switch (event)
            {
            case PHP_EVENT_P_SHORTEN:
                SizeAdjust(true, false);
                return 1;

            case PHP_EVENT_P_RECOVERY:
                if (phragmoplast->DeltaX <= phragmoplast->GetTotalLength())
                {
                    statePlusEnd = PHP_STATE_GROWING; // NEW: move this inside test
                    SizeAdjust(true, true);
                }
                return 1;

            case PHP_EVENT_P_PAUSE:
                statePlusEnd = PHP_STATE_PAUSE;
                return 1;
            }
            break;

        case PHP_STATE_PAUSE:
            switch (event)
            {
            case PHP_EVENT_P_RECOVERY:
                if (phragmoplast->DeltaX <= phragmoplast->GetTotalLength())
                {
                    statePlusEnd = PHP_STATE_GROWING; // NEW: move this inside test
                    SizeAdjust(true, true);
                }
                return 1;

            case PHP_EVENT_P_CATASTROPHE:
                statePlusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(true, false);
                return 1;
            }
            break;

        case PHP_STATE_MATURE:
            if (event == PHP_EVENT_P_CATASTROPHE)
            {
                statePlusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(true, false);
            }
            return 1;

        case PHP_STATE_IN_CP:
            if (event == PHP_EVENT_P_CATASTROPHE)
            {
                statePlusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(true, false);
            }
            return 1;
        case PHP_STATE_IN_DZ:
            if (event == PHP_EVENT_P_CATASTROPHE)
            {
                statePlusEnd = PHP_STATE_SHRINKING;
                SizeAdjust(true, false);
            }

            return 1;
        }
    }
    return 0;
}

/**
 * @brief Selects which event to do next (Monte Carlo) for this MT
 *
 * @return double   time to do next event
 */
double Microtubule::SelectNextEvent()
{
    double dt;
    int event_min = 0;
    double dtt, dt_min = 1e64, L_MT, xp, xm;

    // Polymerisation of + end
    dt_min = 1e32;
    event_min = PHP_EVENT_NONE;

    L_MT = Length(); // length of this MT
    if (towardsDistalZone)
    {
        xp = x2;
        xm = x1;
    }
    else
    {
        xm = x2;
        xp = x1;
    }

    // Polymerisation of - end
    if (typeMT == PHP_TYPE_TREADMILL)
    {
        switch (stateMinusEnd)
        {
        case PHP_STATE_GROWING:
            EVENT(phragmoplast->R_me_polym() / phragmoplast->DeltaX, PHP_EVENT_M_GROW);
            EVENT(phragmoplast->R_me_gs(), PHP_EVENT_M_CATASTROPHE);
            EVENT(phragmoplast->R_me_gp(), PHP_EVENT_M_PAUSE);
            break;
        case PHP_STATE_SHRINKING:
            EVENT(phragmoplast->get_r_me_depolym() / phragmoplast->DeltaX, PHP_EVENT_M_SHORTEN);
            EVENT(phragmoplast->R_me_sg(), PHP_EVENT_M_RECOVERY);
            EVENT(phragmoplast->R_me_sp(), PHP_EVENT_M_PAUSE);
            break;

        case PHP_STATE_PAUSE:
            EVENT(phragmoplast->R_me_pg(), PHP_EVENT_M_RECOVERY);
            EVENT(phragmoplast->R_me_ps(), PHP_EVENT_M_CATASTROPHE);
            break;

        case PHP_STATE_MATURE:
            EVENT(phragmoplast->R_me_gs(), PHP_EVENT_M_CATASTROPHE);
            break;
        default:
            cout << "invalid state : " << statePlusEnd << "\n";
            exit(0);
        }
        dt = dt_min;
    }

    // + end
    switch (statePlusEnd)
    {
    case PHP_STATE_EMPTY: // to improve
        EVENT(phragmoplast->RateReseed(xInit, typeMT), PHP_EVENT_P_RESEED);
        break;

    case PHP_STATE_GROWING:
        EVENT(phragmoplast->RatePolym() / phragmoplast->DeltaX, PHP_EVENT_P_GROW); /// phragmoplast->deltaX
        EVENT(phragmoplast->RateGrowth2Shrink(), PHP_EVENT_P_CATASTROPHE);
        EVENT(phragmoplast->RateGrowth2Pause(), PHP_EVENT_P_PAUSE);
        break;
    case PHP_STATE_SHRINKING:
        EVENT(phragmoplast->GetRateDepolym() / phragmoplast->DeltaX, PHP_EVENT_P_SHORTEN); /// phragmoplast->deltaX
        EVENT(phragmoplast->RateShrink2Growth(), PHP_EVENT_P_RECOVERY);
        EVENT(phragmoplast->RateShrink2Pause(), PHP_EVENT_P_PAUSE);
        break;

    case PHP_STATE_PAUSE:
        EVENT(phragmoplast->RatePause2Growth(), PHP_EVENT_P_RECOVERY);
        EVENT(phragmoplast->RatePause2Shrink(), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_MATURE:
        EVENT(phragmoplast->RateGrowth2Shrink(), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_IN_CP:
        EVENT(phragmoplast->RatePause2ShrinkCellPlate(), PHP_EVENT_P_CATASTROPHE);
        break;

    case PHP_STATE_IN_DZ:
        EVENT(phragmoplast->RatePause2ShrinkDistalZone(), PHP_EVENT_P_CATASTROPHE);
        break;

    default:
        cout << "invalid state : " << statePlusEnd << "\n";
        exit(0);
    }
    dt = dt_min;
    minEvent = event_min;
    return dt;
}

/**
 * @brief Evolve this MT until t_max
 *
 * @param t_max     time to evolve until
 */
void Microtubule::RunUntil(double t_max)
{
    double dt;
    int event;

    while (time < t_max)
    {
        dt = SelectNextEvent();
        DoEvent(minEvent);
        time += dt;
    }
}

/**
 * @brief Bleach this microtubule
 *
 */
void Microtubule::Bleach()
{
    double x, y;
    int i = 0;

    if (typeMT == PHP_TYPE_GROW_INSIDE)
    {
        // swap x1 and x2 temporaily
        x = x1;
        x1 = x2;
        x2 = x;
        y = y1;
        y1 = y2;
        y2 = y;
    }

    if (IsInsideFRAP(x1, y1))
    {
        bleachedx1 = x1;
        ++i;
    }
    else
    {
        // cout<<"x1,y1 Not Inside\n";
    }

    if (IsInsideFRAP(x2, y2))
    {
        if (i == 0)
        {
            bleachedx1 = x2;
        }
        else
        {
            bleachedx2 = x2;
        }
        ++i;
    }
    else
    {
        // cout<<"x2,y2 Not Inside\n";
    }

    if (i < 2)
    {
        if ((y = CrossLeftFRAP()) >= 0)
        {
            if (i == 0)
            {
                bleachedx1 = phragmoplast->FRAPX1;
            }
            else
            {
                bleachedx2 = phragmoplast->FRAPX1;
            }
            ++i;
        }

        if ((y = CrossRightFRAP()) >= 0)
        {
            if (i == 0)
            {
                bleachedx1 = phragmoplast->FRAPX2;
            }
            else
            {
                bleachedx2 = phragmoplast->FRAPX2;
            }
            ++i;
        }

        if ((x = CrossBottomFRAP()) >= 0)
        {
            if (i == 0)
            {
                bleachedx1 = x;
            }
            else
            {
                bleachedx2 = x;
            }
            ++i;
        }

        if ((x = CrossTopFrap()) >= 0)
        {
            if (i == 0)
            {
                bleachedx1 = x;
            }
            else
            {
                bleachedx2 = x;
            }
            ++i;
        }
    }

    if (i >= 2)
    {
        if (bleachedx1 > bleachedx2)
        {
            x = bleachedx1;
            bleachedx1 = bleachedx2;
            bleachedx2 = x;
        }
    }
    else
    {
        bleachedx1 = bleachedx2 = -1;
    }

    if (typeMT == PHP_TYPE_GROW_INSIDE)
    {
        // swap x1 and x2  back
        x = x1;
        x1 = x2;
        x2 = x;
        y = y1;
        y1 = y2;
        y2 = y;
    }
}

/**
 * @brief Save this MT as 1 line in txt file
 *
 * @param ofs   output file stream
 */
void Microtubule::TextSave(ofstream &ofs)
{
    ofs << typeMT << " " << statePlusEnd << " " << x1 << " " << y1 << " " << x2 << " " << y2 << " " << bleachedx1 << " "
        << bleachedx2 << "\n";
}

/**
 * @brief Tests whether a point is inside the FRAP region
 *
 * @param x         x coordinate
 * @param y         y coordinate
 * @return true     if inside
 * @return false    if outside
 */
bool Microtubule::IsInsideFRAP(double x, double y)
{
    return ((phragmoplast->FRAPX1 <= x) && (x <= phragmoplast->FRAPX2) && (phragmoplast->FRAPY1 <= y) &&
            (y <= phragmoplast->FRAPY2));
}

/**
 * @brief Y coordinate of cross with left edge
 *
 * @return double   y coordinate or -1 if none
 */
double Microtubule::CrossLeftFRAP()
{
    double yl = LineY(phragmoplast->FRAPX1);
    if ((((y1 <= yl + 1e-8) && (yl <= y2 + 1e-8)) || ((y2 <= yl + 1e-8) && (yl <= y1 + 1e-8))) &&
        (phragmoplast->FRAPY1 <= yl) && (yl <= phragmoplast->FRAPY2))
    {
        return (yl);
    }
    return (-1);
}

/**
 * @brief Y coordinate of cross with right edge
 *
 * @return double   y coordinate or -1 if none
 */
double Microtubule::CrossRightFRAP()
{
    double yr = LineY(phragmoplast->FRAPX2);
    // std::cout<<" yr="<<yr<<"\n";
    if ((((y1 <= yr + 1e-8) && (yr <= y2 + 1e-8)) || ((y2 <= yr + 1e-8) && (yr <= y1 + 1e-8))) &&
        (phragmoplast->FRAPY1 <= yr) && (yr <= phragmoplast->FRAPY2))
    {
        return (yr);
    }
    return (-1);
}

/**
 * @brief Y coordinate of cross with bottom edge
 *
 * @return double   y coordinate or -1 if none
 */
double Microtubule::CrossBottomFRAP()
{
    double xb = LineX(phragmoplast->FRAPY1);
    // std::cout<<" xb="<<xb<<"\n";
    if ((((x1 <= xb + 1e-8) && (xb <= x2 + 1e-8)) || ((x2 <= xb + 1e-8) && (xb <= x1 + 1e-8))) &&
        (phragmoplast->FRAPX1 <= xb) && (xb <= phragmoplast->FRAPX2))
    {
        return (xb);
    }
    return (-1);
}

/**
 * @brief Y coordinate of cross with top edge
 *
 * @return double   y coordinate or -1 if none
 */
double Microtubule::CrossTopFrap()
{
    double xt = LineX(phragmoplast->FRAPY2);
    // std::cout<<" xt="<<xt<<"\n";
    if ((((x1 <= xt + 1e-8) && (xt <= x2 + 1e-8)) || ((x2 <= xt + 1e-8) && (xt <= x1 + 1e-8))) &&
        (phragmoplast->FRAPX1 <= xt) && (xt <= phragmoplast->FRAPX2))
    {
        return (xt);
    }
    return (-1);
}