#ifndef AMBER_CHRBASEREGION_H
#define AMBER_CHRBASEREGION_H

#include <string>

namespace amber {

// 對應 hmf-common ChrBaseRegion。start/end 皆為 1-based 且含端點。
struct ChrBaseRegion
{
    std::string chromosome;
    int start;
    int end;

    bool containsPosition(const std::string &chr, int position) const {
        return chromosome == chr && position >= start && position <= end;
    }
};

}

#endif
