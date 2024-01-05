#ifndef T_READPARS_HPP
#define T_READPARS_HPP

#include "php_const.hpp"

using namespace std;

enum eParType
{
    T_NONE = 0,
    T_INT = 1,
    T_DOUBLE = 3,
    T_STRING = 4,
};

struct sParameter
{
    const char *name;
    eParType type;
    void *p;
    const char *help;
};



bool ReadOneLine(std::ifstream *ifs, std::string &name, std::string &val, int &line_no);
void ReadPars(const char *fname, struct sParameter *par);
int StringIndex(const char **ptr, const char *name);

#endif
