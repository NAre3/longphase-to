#ifndef AMBER_AMBERSITE_H
#define AMBER_AMBERSITE_H

#include <string>

namespace amber {

struct AmberSite
{
    std::string chromosome;
    int position;
    std::string ref;
    std::string alt;
    bool snpCheck;
    double gnomadFrequency;
};

}

#endif
