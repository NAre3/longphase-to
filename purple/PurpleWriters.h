#ifndef PURPLE_WRITERS_H
#define PURPLE_WRITERS_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

// cnv.somatic.tsv／purity.tsv 寫出的字串（"%.4f"），以及 minorAlleleCopyNumber／majorAlleleCopyNumber 欄的值。
// germline 後置過濾用它們取得與檔案逐位相同的數值。
std::string formatNumber(double value);
double minorAlleleCopyNumber(const PurpleCopyNumber &copyNumber);
double majorAlleleCopyNumber(const PurpleCopyNumber &copyNumber);

void writeCoreOutputs(const std::string &outputDir, const InputData &inputs,
        const std::vector<FittedPurity> &fits, const BestFit &bestFit,
        const std::vector<ObservedRegion> &fittedRegions,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary);

}

#endif
