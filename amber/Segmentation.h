#ifndef AMBER_SEGMENTATION_H
#define AMBER_SEGMENTATION_H

#include <string>
#include <vector>

#include "AmberOutput.h"

namespace amber {

// 對應 hmf-common 的 common.segmentation 與 common.segmentation.copynumber 套件，
// 以及 amber 的 BAFSegmenter。行為出處逐條見 Segmentation.cpp 的註解。

// BAFSegmenter 建構時傳入 gamma = 100.0（BAFSegmenter.java:37），AMBER 4.3 無 CLI 可調
constexpr double BAF_SEGMENTATION_GAMMA = 100.0;

// PerArmSegmenter 的 uniform-penalty 門檻。AMBER 4.3 呼叫端使用的值（見 CP-A8 的
// uniformPenaltyThreshold 欄）。本研究的 dev 樣本 totalCount = 701544 遠高於此，
// 走 per-arm-gamma 分支；uniform 分支未被驗證（spec U5）。
constexpr int UNIFORM_PENALTY_THRESHOLD = 100000;

struct PcfSegment
{
    std::string chromosome; // shortName，例如 "1"（不含 chr 前綴）
    int start = 0;
    int end = 0;
    double meanRatio = 0.0;
};

struct ArmSegments
{
    std::string armId;              // 例如 "1_P"
    std::string chromosomeShort;    // 例如 "1"
    int chromosomeRank = 0;         // HumanChromosome 的 enum 序，用於排序
    char arm = 'P';
    std::vector<double> values;     // valuesForSegmentation（= rawValues，BAFSegmenter 兩者相同）
    std::vector<PcfSegment> segments;
};

struct SegmentationResult
{
    int totalCount = 0;
    int uniformPenaltyThreshold = UNIFORM_PENALTY_THRESHOLD;
    std::string penaltyMode;  // "uniform" 或 "per-arm-gamma"
    double gamma = BAF_SEGMENTATION_GAMMA;
    std::vector<ArmSegments> arms; // 已依 ChrArm.compareTo 排序
};

// 對應 BAFSegmenter.writeSegments → PerArmSegmenter 建構 + getSegmentation。
// value(baf) = tumorModifiedBAF，入選條件為 value >= 0.0。
SegmentationResult segmentBafs(const std::vector<AmberBAF> &bafs);

// 對應 SegmentsFile.write（chromosome 加 chr 前綴、MeanRatio 以 DecimalFormat("#.####") 格式化）
void writeSegmentsFile(const std::string &path, const SegmentationResult &result);

// CP-A8 / CP-A9 的 dump。與 Java 端 patch_segmenter.py 插入的位置對應。
void writeSegmentationCheckpoints(const SegmentationResult &result);

// 對外暴露供測試：Runmed（含 smoothEnds）與 Gamma 的 segment penalty
std::vector<double> runmed(const std::vector<double> &data, int k, bool smooth);
double gammaSegmentPenalty(const std::vector<double> &y, double gamma, bool normalise);

}

#endif
