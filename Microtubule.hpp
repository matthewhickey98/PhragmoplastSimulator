#ifndef MICROTUBLE_HPP
#define MICROTUBLE_HPP

#include "Bitmap.hpp"
#include "Phragmoplast.hpp"
#include "php_const.hpp"

class Phragmoplast;

// A single microtubule

/**
 * @brief A class representing a single microtubule.
 * 
 */
class Microtubule
{
  protected:
    double xInit;               // Initial X location
    double yInit;               // Initial Y location
    double x1;                  // Position of the seed
    double x2;                  // End of growing MT
    double y1;                  // Position of the seed
    double y2;                  // End of growing MT
    double deltaX;              // X increment
    double deltaY;              // Y increment
    double theta;               // Orientation
    double bleachedx1;          // Left coordinate bleached region
    double bleachedx2;          // Right coordinate bleached region
    double time;                // time
    int typeMT;                 // TREADMILL, SEEDED (GROWING_OUTSIDE, GROWING_INSIDE)
    bool towardsDistalZone;     // true -> +end away from CP
    int statePlusEnd;           // +end : EMPTY, GROWING, SHRINKING, MATURE
    int stateMinusEnd;          // -end : EMPTY, GROWING, SHRINKING, MATURE
    int minEvent;               // Event for this microtubule
    Phragmoplast *phragmoplast; // Pointer to phragmoplast

  public:
    Microtubule() : bleachedx1(-1), bleachedx2(-1), typeMT(PHP_TYPE_NONE), minEvent(0), phragmoplast(0) {};
    ~Microtubule(){

    };

    bool IsInsideFRAP(double x, double y);

    void Init(double x1, double y1, double L, double theta, int type, Phragmoplast *php);
    void Draw(Bitmap &bm);
    void Bleach();
    void TextSave(std::ofstream &ofs);
    void CheckMTPosition();
    void UpdateBleach();
    void SizeAdjust(bool plus = true, bool grow = true);
    void RunUntil(double t_max);

    double CrossLeftFRAP();
    double CrossRightFRAP();
    double CrossBottomFRAP();
    double CrossTopFrap();

    int DoEvent(int event);
    double SelectNextEvent();

    inline int GetMinEvent() { return minEvent; };

    /**
     * @brief Get the length of the microtubule.
     * 
     * @return double       the length
     */
    inline double Length()
    {
        return ((this->statePlusEnd == PHP_STATE_EMPTY) ? 0 : sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)));
    };

    // corrdinate of the MT
    inline double LineX(double y)
    {
        return (((x2 - x1) * y + x1 * y2 - x2 * y1) / (y2 - y1));
    }
    inline double LineY(double x)
    {
        return (((y2 - y1) * x + y1 * x2 - y2 * x1) / (x2 - x1));
    }
    inline double GetInitX()
    {
        return xInit;
    };
    inline double GetInitY()
    {
        return yInit;
    };
    inline double GetBleachedX1()
    {
        return bleachedx1;
    };
    inline double GetBleachedX2()
    {
        return bleachedx2;
    };
    inline double GetX1()
    {
        return (x1);
    };
    inline double GetX2()
    {
        return (x2);
    };
    inline double GetY1()
    {
        return (y1);
    };
    inline double GetY2()
    {
        return (y2);
    };

    inline double GetTime()
    {
        return (time);
    }
    inline int GetType()
    {
        return (typeMT);
    }
    inline int GetPlusState()
    {
        return (statePlusEnd);
    }
    inline int GetMinusState()
    {
        return (stateMinusEnd);
    }

    inline Microtubule &SetTime(double t)
    {
        time = t;
        return *this;
    }
    inline void SetType(int type)
    {
        typeMT = type;
    }
    inline void SetPlusState(int state)
    {
        statePlusEnd = state;
    }
    inline void SetMinusState(int state)
    {
        stateMinusEnd = state;
    }
};
#endif
