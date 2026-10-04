#ifndef PURPLE_OBSERVED_H
#define PURPLE_OBSERVED_H

#include <vector>

#include "PurpleTypes.h"

namespace purple {

std::vector<ObservedRegion> createObservedRegions(const InputData &data, const std::vector<SupportSegment> &segments);
void dumpObservedRegions(const std::vector<ObservedRegion> &regions);

}

#endif
