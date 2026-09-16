#ifndef PURPLE_SUMMARY_H
#define PURPLE_SUMMARY_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

SummaryContext buildSummaryContext(const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const std::string &ensemblDataDir);
void dumpSummaryContext(const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary);

}

#endif
