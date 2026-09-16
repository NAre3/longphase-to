#ifndef PURPLE_COPY_NUMBER_H
#define PURPLE_COPY_NUMBER_H
#include <vector>
#include "PurpleTypes.h"
namespace purple {
std::vector<PurpleCopyNumber> buildCopyNumbers(const std::vector<ObservedRegion> &fittedRegions, const FittedPurity &fit);
void dumpCopyNumbers(const std::vector<PurpleCopyNumber> &copyNumbers);
}
#endif
