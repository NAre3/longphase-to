#ifndef PURPLE_FITTING_H
#define PURPLE_FITTING_H

#include <vector>
#include "PurpleTypes.h"

namespace purple {

std::vector<const ObservedRegion *> selectFittingRegions(const std::vector<ObservedRegion> &regions);
std::vector<FittedPurity> fitPurityGrid(const std::vector<const ObservedRegion *> &regions, int averageTumorDepth, int threads);
void dumpFittingRegions(const std::vector<const ObservedRegion *> &regions);
void dumpPurityGrid(const std::vector<FittedPurity> &fits);

}

#endif
