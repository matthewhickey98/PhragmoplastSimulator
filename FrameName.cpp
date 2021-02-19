#include "FrameName.h"
#include <string.h>
#include <stdio.h>
#include <iostream>

FrameName::FrameName(const char *prefix, const char *suffix) : index_(-1)
{
    strncpy(prefix_, prefix, MAX_FILE_SIZE - 1);
    strncpy(suffix_, suffix, MAX_FILE_SIZE - 1);
    //std::cerr<<"CONST: index_="<<index_<<"\n";
};

char *FrameName::next_name()
{
    inc();
    //std::cerr<<">>> index_="<<index_<<"\n";
    sprintf(buffer_, "%s_%4.4d%s", prefix_, index_, suffix_);
    return (buffer_);
};

char *FrameName::current_name_other_suffix(const char *suffix)
{
    sprintf(buffer2_, "%s_%4.4d%s", prefix_, index_, suffix);
    return (buffer2_);
};

void FrameName::set_prefix(const char *prefix)
{
    strncpy(prefix_, prefix, MAX_FILE_SIZE - 1);
}
void FrameName::set_suffix(const char *suffix)
{
    strncpy(suffix_, suffix, MAX_FILE_SIZE - 1);
}
