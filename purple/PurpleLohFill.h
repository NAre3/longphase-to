#ifndef PURPLE_PURPLELOHFILL_H
#define PURPLE_PURPLELOHFILL_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

// ---- 用 longphase-to 的 LOH 補回 AMBER 在高純度漏掉的 BAF 觀測 ----
//
// 純度接近 1.0 時，真正的 somatic LOH 位點只剩一個 allele，AMBER 無法確認它原本是
// germline 雜合而把它篩掉，LOH 段的 bafCount 幾乎歸零、殘留位點的 BAF 也偏低，
// PURPLE 因此失去排除低純度解的證據。這裡在 AMBER 自己也缺資料的段上，
// 把 bafCount 補成「AMBER 沒篩掉的話應有的位點數」、observedBAF 設為 1.0。
//
// 規則：
//   ok   = 體染色體 && DIPLOID && 0 <= observedTumorRatio <= 3（不要求 bafCount > 0）
//   cov  = 段與 LOH 區間的重疊長度 / 段長（兩者都當閉區間）
//   pan  = raw germline site panel 落在 [start, end] 的位置數（不套 tumor-only blacklist）
//   ret  = Σ bafCount / Σ pan，取 ok && cov < 0.1 && bafCount > 0 的段
//   若 ok && cov >= 0.9 && pan > 0 && bafCount / (ret * pan) < 0.3
//      則 bafCount = max(1, round-half-even(ret * pan))，observedBAF = 1.0
//
// key 一律用去掉 chr 前綴後的染色體名。
struct LohFillInput
{
    std::map<std::string, std::vector<std::pair<int, int>>> lohSpans;
    std::map<std::string, std::vector<int>> panelPositions;  // 已排序
};

struct LohFillRecord
{
    std::string chromosome;
    int start = 0;
    int end = 0;
    int originalBafCount = 0;
    double originalObservedBaf = 0;
    int panelCount = 0;
    double coverage = 0;
    double expected = 0;   // E2 = ret * pan
    double retention = 0;  // r = bafCount / E2
    int filledBafCount = 0;
};

struct LohFillResult
{
    double ret = 0;
    int retSegments = 0;
    std::vector<LohFillRecord> filled;
};

// 只收體染色體（1–22）。
std::map<std::string, std::vector<int>> loadPanelPositions(const std::string &path);
std::map<std::string, std::vector<std::pair<int, int>>> loadLohBed(const std::string &path);

LohFillResult applyLohFill(std::vector<ObservedRegion> &regions, const LohFillInput &input);

// <sample>.purple.loh_fill.tsv：被補段的原始值與補回值，供追溯。
void writeLohFillTsv(const std::string &outputDir, const std::string &sampleId, const LohFillResult &result);
void dumpLohFill(const LohFillResult &result);

}

#endif
