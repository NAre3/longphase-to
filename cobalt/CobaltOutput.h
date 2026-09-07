#ifndef COBALT_COBALTOUTPUT_H
#define COBALT_COBALTOUTPUT_H

#include <string>
#include <vector>

#include "BamRatio.h"
#include "GcBuckets.h"

namespace cobalt {

// CobaltRatio.java 的 tumor-only 版本（四個 reference 欄位固定 -1.0，G14）
struct CobaltRatio
{
    std::string chromosome;      // 帶 chr 前綴（versionedChromosome）
    int position = 0;
    double referenceReadDepth = -1.0;
    double referenceGCRatio = -1.0;
    double referenceGCContent = -1.0;
    double referenceGCDiploidRatio = -1.0;
    double tumorReadDepth = 0.0;
    double tumorGCRatio = 0.0;
    double tumorGcContent = 0.0;
};

std::vector<CobaltRatio> collateResults(const std::vector<BamRatio> &tumorRatios);

// CobaltRatioFile.write：gzip、DecimalFormat("#.####")、**檔案欄序 != 物件欄序**
void writeCobaltRatioFile(const std::string &path, const std::vector<CobaltRatio> &ratios);

// CobaltGcMedianFile.write：%.2f，只輸出 MeanDepths[i] > 0 的 bucket
void writeGcMedianFile(const std::string &path, double sampleMean, double sampleMedian,
                       const GcBucketStatistics &stats);

}

#endif
