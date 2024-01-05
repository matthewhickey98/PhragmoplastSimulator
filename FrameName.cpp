#include "FrameName.hpp"

/**
 * @brief Construct a new Frame Name
 * 
 * @param pre   Prefix of the frame name
 * @param suf   Suffix of the frame name
 */
FrameName::FrameName(const std::string pre, const std::string suf) : index(-1)
{
    prefix = pre;
    suffix = suf;
};

/**
 * @brief Get the next name
 * 
 * @return const std::string&   Next name
 */
const std::string &FrameName::NextName()
{
    Increment();
    std::stringstream ss;
    ss << prefix;
    ss << index;
    ss << suffix;
    buffer = ss.str();
    return buffer;
};

/**
 * @brief Get the Current Name with a different suffix
 * 
 * @param suf   Suffix of the frame name
 * @return      const std::string&   Name
 */
const std::string &FrameName::CurrentNameOtherSuffix(const std::string &suf)
{
    std::stringstream ss;
    ss << prefix;
    ss << index;
    ss << suf;
    buffer2 = ss.str();
    return (buffer2);
};

