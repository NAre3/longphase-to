#ifndef COBALT_SEGMENTATION_H
#define COBALT_SEGMENTATION_H

#include <string>
#include <vector>

#include "CobaltOutput.h"
#include "../common/Segmentation.h"

namespace cobalt {

// PerArmSegmenter 的 uniform 分支門檻（PerArmSegmenter.java:39）
constexpr int UNIFORM_PENALTY_THRESHOLD = 100000;

struct PcfSegmentOut
{
    std::string chromosome;      // 短名
    int start = 0;
    int end = 0;
    double meanRatio = 0.0;      // rawValues 該段的平均（**非** pcfMeans）
};

// DataForSegmentation：**兩個陣列**（COBALT 的 valuesForSegmentation 與 rawValues 不同，
// 這正是 G24f 指出 AMBER 的單一 values 陣列無法表達的地方）
struct ArmData
{
    std::string armId;             // 例如 "1_P"（dump 用 "ChrArm[chromosome=1, arm=P]"）
    std::string chromosomeShort;
    int chromosomeRank = 0;
    char arm = 'P';
    std::vector<double> valuesForSegmentation;   // (float) 量化後的 log2
    std::vector<double> rawValues;               // 未經 floor 的原 ratio
    std::vector<int> positions;                  // 各值對應的 window position

    // 中間量與結果
    lp::GammaTrace trace;
    lp::Fit fit;
    std::vector<double> pcfMeans;                // PiecewiseConstantFit 自己的 means[]
    std::vector<PcfSegmentOut> segments;
};

struct SegmentationResult
{
    int totalCount = 0;
    int uniformPenaltyThreshold = UNIFORM_PENALTY_THRESHOLD;
    std::string penaltyMode;     // "uniform" 或 "per-arm-gamma"
    double gamma = 100.0;
    bool isWindowed = true;
    std::vector<ArmData> arms;   // 依 ChrArm.compareTo 排序
};

SegmentationResult segmentRatios(const std::vector<CobaltRatio> &ratios, double gamma);

// SegmentsFile.write：表頭四欄 Chromosome/Start/End/MeanRatio（無 n.probes，findings F1）
void writeSegmentsFile(const std::string &path, const SegmentationResult &result);

}

#endif
