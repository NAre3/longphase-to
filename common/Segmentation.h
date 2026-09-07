#ifndef LP_SEGMENTATION_H
#define LP_SEGMENTATION_H

#include <string>
#include <vector>

namespace lp {

// hmf-common 的 common.segmentation.copynumber 通用核心：Runmed、Gamma penalty 與
// Segmenter 的最小成本分段。**這一層不知道 BAF 或 ratio**——AMBER 與 COBALT 各自
// 提供 value(...) 與輸出格式。行為出處逐條見 Segmentation.cpp 的註解。

// ChrArmLocator：以 position < centromere 判定 P 臂（GRCh38 座標表見 Centromeres38.inc）
int centromereOf(const std::string &chromosome);

// HumanChromosome 的宣告順序（_1.._22, _X, _Y），ChrArm.compareTo 的第一鍵
int chromosomeOrdinal(const std::string &shortName);

// Doubles.mean：依序累加後除以長度。累加順序影響最後一位，不要改成 pairwise 或 Kahan。
double mean(const double *values, std::size_t count);

// Runmed（Runmed.java）：k 須為奇數；smooth 對應 smoothEnds
std::vector<double> runmed(const std::vector<double> &data, int k, bool smooth);

// Gamma（Gamma.java）：normalise 為假時直接回傳 gamma
double gammaSegmentPenalty(const std::vector<double> &y, double gamma, bool normalise);

// 【EXP-C008 新增】gammaSegmentPenalty 的中間量，供 COBALT 的 CP-C13b dump。
// **既有的三參數版本原樣保留、行為不變**（轉呼叫此版並傳 nullptr），
// 故 AMBER 的呼叫點一行不改、執行路徑逐位元相同。
struct GammaTrace
{
    int n = 0;
    int filterWidth = 0;
    std::vector<double> runningMedians;
    double mad = 0.0;
    bool normalised = false;
    double penalty = 0.0;
};

double gammaSegmentPenalty(const std::vector<double> &y, double gamma, bool normalise, GammaTrace *trace);

// PiecewiseConstantFit 的分段結果：各段長度與起始索引
struct Fit
{
    std::vector<int> lengths;
    std::vector<int> startPositions;
};

// Segmenter（Segmenter.java）：最小成本分段的動態規劃
Fit segment(const std::vector<double> &y, double segmentPenalty);

// 【EXP-C008 新增】PiecewiseConstantFit 自己的 means[]，供 COBALT 的 CP-C13c dump。
// means[i] = round(Doubles.mean(segment))，round(v) = Math.round(v*1000)/1000.0
// （Math.round(double) 規格為 (long) floor(x + 0.5)，C++ 端須寫 floor(x+0.5)，不可用 llround）
// （Segmentation.java:48,142-145）。該值隨後被 rawValues 的平均取代，只有在此觀察得到。
// **既有的兩參數版本原樣保留、行為不變。**
Fit segment(const std::vector<double> &y, double segmentPenalty, std::vector<double> *pcfMeans);

}

#endif
