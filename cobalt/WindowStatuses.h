#ifndef COBALT_WINDOWSTATUSES_H
#define COBALT_WINDOWSTATUSES_H

#include <string>
#include <vector>

#include "GcProfile.h"
#include "Regions.h"

namespace cobalt {

// WindowStatus.java
struct WindowStatus
{
    std::string chromosome;      // gcProfile.chromosome()，帶 chr 前綴
    int start = 0;
    int end = 0;
    bool excluded = false;
    bool unmappable = false;
    bool nonDiploid = false;

    bool maskedOut() const { return excluded || unmappable || nonDiploid; }
};

struct WindowStatuses
{
    std::vector<std::string> chromosomes;                 // 短名，依 ordinal 昇冪
    std::vector<std::vector<WindowStatus>> byChromosome;
};

WindowStatuses buildWindowStatuses(const GcProfileData &gcData,
                                   const std::vector<ExcludedRegion> &exclusions,
                                   const DiploidRegions &diploidRegions);

}

#endif
