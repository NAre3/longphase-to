#ifndef PURPLE_PURPLEPHASEDBAF_H
#define PURPLE_PURPLEPHASEDBAF_H

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

// ---- 以 phased SNV 的 haplotype 方向 HMM 估擬合段的 BAF ----
//
// 原版段 BAF 是逐位點折疊 0.5 + |BAF − 0.5| 後取中位數。低純度時真實 BAF 貼近
// 0.5、單點雜訊又大，折疊會把雜訊翻成往上的系統偏差（平衡段約 +0.04～+0.07）。
// 這裡改用 longphase-to 自己的 phasing：位點換算成 hap1 的 read 比例
// h = (hap1 帶 alt) ? tumorBaf : 1 − tumorBaf，先依 haplotype 合併、最後才折疊一次。
//
// phase block 內仍有方向翻轉（switch error 或段內換邊），整個 block 直接合併會把
// 不平衡段拉向 0.5，逐格依正負翻轉又在低純度把折疊偏差帶回來。所以把方向當成
// 隱藏狀態，對方向邊際化：
//   小格 c ＝ (PS, 位置 / 100 kb)；H_c ＝ Σ d·h / Σ d，n_c ＝ Σ d；n_c < 20 的小格丟棄；
//     小格依第一個位點的位置排序。
//   狀態 +（hap1 為多數）／−；發射 N(f 或 1 − f, f(1 − f)·κ / n_c)，κ 修正長讀
//     一條 read 被多個位點重複計入的變異膨脹。
//   轉移：相鄰兩格同一 PS 時以機率 s 翻轉；PS 不同時方向重新各 0.5。
//   段 BAF ＝ f 在 0.5～0.995（間距 0.0025）上的最大概似值（前向演算法）。
//   κ ∈ {1,2,4,8,12,16,24,32}、s ∈ {0.001,0.003,0.01,0.03,0.1}：取所有段
//     max_f 對數概似總和最大者（κ 外圈、s 內圈，相同時取先出現者）。
// 沒有可用小格的段保留原 observedBAF。LOH 補回的段保留補回值（仍參與 κ、s 的選擇）。
//
// 已知相近方法：Numbat／hahmmr 的 haplotype-aware allele HMM。

// key：去掉 chr 前綴的染色體 → 1-based 位置 → (hap1 帶 alt, PS)
using PhasedSnvMap = std::map<std::string, std::unordered_map<int, std::pair<bool, int>>>;

struct PhasedBafResult
{
    double kappa = 0;
    double switchProbability = 0;
    int regionsWithCells = 0;
    int regionsUpdated = 0;
    // 與 regions 同序：原 observedBaf、HMM 值（未更新則為原值）、可用小格數、是否為 LOH 補回段
    std::vector<double> originalBaf;
    std::vector<double> phasedBaf;
    std::vector<int> cellCount;
    std::vector<char> lohFilled;
};

// regions：要更新的擬合段（呼叫端的副本，observedBaf 會被改寫）。
// lohFilledKeys：LOH 補回段的 (去 chr 前綴的染色體, start)。
PhasedBafResult applyPhasedBaf(std::vector<ObservedRegion *> &regions, const std::vector<AmberBaf> &bafs,
        const PhasedSnvMap &phased, const std::set<std::pair<std::string, int>> &lohFilledKeys);

void dumpPhasedBaf(const std::vector<ObservedRegion *> &regions, const PhasedBafResult &result);

}

#endif
