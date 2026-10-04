#ifndef PURPLE_FITTING_H
#define PURPLE_FITTING_H

#include <vector>
#include "PurpleTypes.h"

namespace purple {

std::vector<const ObservedRegion *> selectFittingRegions(const std::vector<ObservedRegion> &regions);
std::vector<FittedPurity> fitPurityGrid(const std::vector<const ObservedRegion *> &regions, int averageTumorDepth, int threads);
BestFit selectTumorOnlyBestFit(const std::vector<FittedPurity> &fits, const std::vector<ObservedRegion> &regions);
std::vector<ObservedRegion> fitObservedRegions(
        const std::vector<ObservedRegion> &regions, const FittedPurity &fit, int averageTumorDepth, Gender gender);
void dumpFittingRegions(const std::vector<const ObservedRegion *> &regions);
void dumpPurityGrid(const std::vector<FittedPurity> &fits);
void dumpBestFit(const BestFit &best);
void dumpFittedRegions(const std::vector<ObservedRegion> &regions);

}

#endif
