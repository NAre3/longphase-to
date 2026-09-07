#ifndef COBALT_CONSOLIDATION_H
#define COBALT_CONSOLIDATION_H

#include <string>
#include <vector>

#include "BamRatio.h"

namespace cobalt {

// LowCovBucket.java
struct LowCovBucket
{
    int startPosition = 0;
    int endPosition = 0;
    int bucketPosition = 0;
};

// ResultsConsolidator.calcConsolidationCount
int calcConsolidationCount(double medianReadDepth);

// ResultsConsolidator.consolidateIntoBuckets
std::vector<LowCovBucket> consolidateIntoBuckets(const std::vector<int> &windowPositions, int consolidationCount);

// 對應 WholeGenome.resultsConsolidator 的分支選擇結果。
// NoOp 時 consolidationCount 記為 1。
struct ConsolidatorChoice
{
    std::string className;      // "NoOpConsolidator" 或 "LowCoverageConsolidator"
    int consolidationCount = 1;
};

ConsolidatorChoice chooseConsolidator(double medianReadDepth);

// LowCoverageConsolidator.consolidate + BamRatios.consolidate 的回填。
// ratios 依 (chromosome ordinal, position) 排序、逐染色體連續。就地更新。
// **本函式在 dev 樣本上不會被呼叫**（NoOp 分支），見 behaviour-contract.md §1.1。
void applyLowCoverageConsolidation(std::vector<BamRatio> &ratios, int consolidationCount);

}

#endif
