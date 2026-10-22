#ifndef PURPLE_PURPLELOHTDF_H
#define PURPLE_PURPLELOHTDF_H

#include <string>
#include <vector>

#include "PurplePhasedBaf.h"
#include "PurpleTypes.h"

namespace purple {

// ---- 以 LOH 段估 tumor DNA fraction，限制 PURPLE 的候選範圍 ----
//
// 段的 tumor DNA fraction：f = 1 − 2R(1 − B)，B 為段的主要 haplotype 比例，R 為段深度比 ÷
// chr1–22 擬合段的長度加權平均。LOH 段（minor CN 0）的 f 等於整體 tumor DNA fraction，
// 與 ploidy 無關；其他段為 f·(1 − 2·minor/ploidy)，一定較小（ASCAT 的 nB 式令 nB = 0 即得）。
// 所以 f 分布最右邊的峰是 LOH 群。
//
// 低純度時 PURPLE 會把真 LOH 段解釋成 ploidy 加倍下 minor CN ≈ 1 的段（整數，偏差罰分不罰），
// 選到 ploidy 約兩倍、tumor DNA fraction 偏高的解；正解仍在候選中，分數只差數個百分點。
// 這裡用 LOH 群的 f 只框出範圍，ploidy 與 copy number 仍由 PURPLE 的分數決定：
//   取 chr1–22 擬合段：
//     LOH 補回段（PurpleLohFill.h）→ B = 1.0，權重 = 補回後的 bafCount；
//     其餘有 HMM 小格的段（PurplePhasedBaf.h）→ B = HMM 段 BAF，權重 = bafCount；
//     只留 0 < f ≤ 1。
//   取峰：f 以 0.005 分格、權重加總；每個有資料的格加總左右各 3 格；
//     「峰」＝ 平滑值 ≥ 最高值的 25% 且 ≥ 左右鄰格（鄰格無資料時視為 0），取最右邊的峰，
//     f0 = (格 + 0.5) × 0.005。
//   候選：只保留 |p·ψ / (p·ψ + 2(1 − p)) − f0| ≤ 0.05 者（ψ 為擬合的平均 ploidy），
//     依原分數順序交給原本的選解；範圍內沒有候選時不限制。
// 純腫瘤的 LOH 區沒有 het，HMM 段 BAF 取不到，由 LOH 補回段提供 f = 1 的證據。

struct LohTdfSegment
{
    std::string chromosome;
    int start = 0;
    int end = 0;
    int weight = 0;
    double ratio = 0;   // observedTumorRatio ÷ chr1–22 長度加權平均
    double baf = 0;
    double f = 0;
    bool filled = false;
};

struct LohTdfEstimate
{
    double meanRatio = 0;
    double f0 = 0;
    bool valid = false;
    std::vector<LohTdfSegment> segments;   // 0 < f ≤ 1 的段，依擬合段順序
};

// regions、phased 與 applyPhasedBaf 的輸入輸出同序。
LohTdfEstimate estimateLohTdf(const std::vector<ObservedRegion *> &regions, const PhasedBafResult &phased);

// 依 fits 原順序保留 tumor DNA fraction 落在 f0 ± halfWidth 的候選。
std::vector<FittedPurity> restrictToTdfWindow(const std::vector<FittedPurity> &fits, double f0, double halfWidth);

double tumorDnaFraction(double purity, double ploidy);

void writeLohTdfTsv(const std::string &outputDir, const std::string &sampleId, const LohTdfEstimate &estimate,
        double halfWidth, size_t candidatesKept, size_t candidatesTotal);

}

#endif
