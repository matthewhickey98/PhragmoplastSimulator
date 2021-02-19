#ifndef T_READPARS_H
#define T_READPARS_H
#include <string>
#include <sstream>
#include <fstream>
#include <istream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <string.h>
#include <exception>
using namespace std;
enum ParType
{
    T_NONE = 0,
    T_INT = 1,
    T_DOUBLE = 3,
    T_STRING = 4,
};

struct Parameter
{
    const char *name;
    ParType type;
    void *p;
    const char *help;
};

bool read_one_line(std::ifstream *ifs, std::string &name, std::string &val,
                   int &line_no);
void read_pars(const char *fname, struct Parameter *par);

int string_index(const char **ptr, const char *name);

#endif
