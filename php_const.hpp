#ifndef PHP_CONST_HPP
#define PHP_CONST_HPP

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <exception>
#include <execution>
#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <list>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <random>

const int PHP_STATE_EMPTY = 1;
const int PHP_STATE_GROWING = 2;
const int PHP_STATE_SHRINKING = 3;
const int PHP_STATE_PAUSE = 4;
const int PHP_STATE_MATURE = 5;
const int PHP_STATE_IN_CP = 6; // inside CP
const int PHP_STATE_IN_DZ = 7; // inside DZ

const int PHP_EVENT_MASK_MINUS = 1 << 8;
const int PHP_EVENT_NONE = 0;
const int PHP_EVENT_P_GROW = 1;
const int PHP_EVENT_P_SHORTEN = 2;
const int PHP_EVENT_P_CATASTROPHE = 3;
const int PHP_EVENT_P_PAUSE = 4;
const int PHP_EVENT_P_RECOVERY = 5;
const int PHP_EVENT_P_RESEED = 6;
const int PHP_EVENT_M_GROW = PHP_EVENT_P_GROW | PHP_EVENT_MASK_MINUS;
const int PHP_EVENT_M_SHORTEN = PHP_EVENT_P_SHORTEN | PHP_EVENT_MASK_MINUS;
const int PHP_EVENT_M_CATASTROPHE = PHP_EVENT_P_CATASTROPHE | PHP_EVENT_MASK_MINUS;
const int PHP_EVENT_M_PAUSE = PHP_EVENT_P_PAUSE | PHP_EVENT_MASK_MINUS;
const int PHP_EVENT_M_RECOVERY = PHP_EVENT_P_RECOVERY | PHP_EVENT_MASK_MINUS;
const int PHP_EVENT_M_RESEED = PHP_EVENT_P_RESEED | PHP_EVENT_MASK_MINUS;

const int PHP_TYPE_NONE = 0;
const int PHP_TYPE_TREADMILL = 1;    // + and - end
const int PHP_TYPE_SEEDED = 2;       // + end only
const int PHP_TYPE_GROW_OUTSIDE = 3; // + end -> DZ
const int PHP_TYPE_GROW_INSIDE = 4;  // + end -> PZ

const int PHP_PATH_LENGTH = 512;

#endif
