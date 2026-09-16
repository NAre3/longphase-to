#ifndef PURPLE_WRITERS_H
#define PURPLE_WRITERS_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

void writeCoreOutputs(const std::string &outputDir, const InputData &inputs,
        const std::vector<FittedPurity> &fits, const BestFit &bestFit,
        const std::vector<ObservedRegion> &fittedRegions,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary);

}

#endif
