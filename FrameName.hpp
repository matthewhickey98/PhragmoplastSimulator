#ifndef FNAME_HPP
#define FNAME_HPP
#include "php_const.hpp"

class FrameName
{
public:
    FrameName() : index(-1) {};
    FrameName(const std::string pre, const std::string suf);
    ~FrameName() {};

    const std::string &NextName();
    const std::string &CurrentNameOtherSuffix(const std::string &suf);
    
    inline void SetPrefix(const std::string &pre) { prefix = pre; };
    inline void SetSuffix(const std::string &suf) { suffix = suf; };

    inline const std::string &CurrentName() const { return (buffer); };
    inline void Increment() { ++index; };
    inline void SetIndex(int indexIn) { index = indexIn; };
    inline int GetIndex() { return (index); };

protected:
    std::string prefix;
    std::string suffix;
    std::string buffer;
    std::string buffer2;
    int index;
};

#endif
