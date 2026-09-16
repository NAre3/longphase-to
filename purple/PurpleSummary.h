#ifndef PURPLE_SUMMARY_H
#define PURPLE_SUMMARY_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

void dumpSummaryContext(const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const std::string &ensemblDataDir);

}

#endif
