#ifndef PURPLE_FITTING_H
#define PURPLE_FITTING_H

#include <vector>
#include "PurpleTypes.h"

namespace purple {

// excludeAllosomes：tumor-only 不以 X、Y 擬合（PURPLE 4.0～4.3 的 isAllosome 條件、官方 README 的
// tumor-only 說明）。預設 false＝現行 4.4 行為（只排除 Y）。
std::vector<const ObservedRegion *> selectFittingRegions(const std::vector<ObservedRegion> &regions, bool excludeAllosomes = false);
// ambiguousOverride ≥ 0 時取代 expectedBaf(averageTumorDepth)。phased 段 BAF 只在 haplotype 層級
// 折疊一次，平衡段的期望值就是 0.5（PurplePhasedBaf.h）。
std::vector<FittedPurity> fitPurityGrid(const std::vector<const ObservedRegion *> &regions, int averageTumorDepth, int threads,
        double ambiguousOverride = -1);
BestFit selectTumorOnlyBestFit(const std::vector<FittedPurity> &fits, const std::vector<ObservedRegion> &regions);
// halfAmbiguous 非 nullptr 時與 regions 同序；為 1 的段以 0.5 取代 expectedBaf（phased 段 BAF 的擬合段）。
std::vector<ObservedRegion> fitObservedRegions(
        const std::vector<ObservedRegion> &regions, const FittedPurity &fit, int averageTumorDepth, Gender gender,
        const std::vector<char> *halfAmbiguous = nullptr);
void dumpFittingRegions(const std::vector<const ObservedRegion *> &regions);
void dumpPurityGrid(const std::vector<FittedPurity> &fits);
void dumpBestFit(const BestFit &best);
void dumpFittedRegions(const std::vector<ObservedRegion> &regions);

}

#endif
