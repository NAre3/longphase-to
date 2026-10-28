#ifndef PURPLE_PURPLESCALECHECK_H
#define PURPLE_PURPLESCALECHECK_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

// ---- 段格點尺度檢查與改值 ----
//
// TDF 固定時 ψ 與 2ψ 的 ratio、BAF 期望值相同，PURPLE 的擬合分數只能靠事件懲罰分辨兩者，
// 而事件懲罰偏好較小的尺度，真值為 2ψ 的樣本可能被選成 ψ。這裡看「哪個尺度讓段的
// major／minor 落在整數格點上」：
//   y = 值 / ψ；A(s) = ½(|Σw·e^{2πi·s·y_major}| + |Σw·e^{2πi·s·y_minor}|) / Σw，
//   s ∈ [1.5, 10] 每 0.01，s* = argmax A（相同時取最小 s）。
//   段：自體染色體、長度 ≥ 1 Mb、bafCount ≥ 20；w = 長度。
//   |s*/(2ψ) − 1| ≤ 0.25 → "double"；|s*/ψ − 1| ≤ 0.25 → "keep"；否則 "unclear"。
// 判 "double" 且 ψ ≥ 2 時，改選完整格點中 |TDF − f0| ≤ 0.05、ploidy ∈ [1.75ψ, 2.25ψ] 的最低分候選。
// ψ ≥ 2 是為了避開 hypodiploid（ψ < 2）且格點乾淨時 A(ψ) ≈ A(2ψ) 的誤判。
//
// ψ 與 major／minor 先經過與 purity.tsv、cnv.somatic.tsv 相同的 "%.4f" 捨入再計算，
// 讓結果與「從輸出檔重算」一致。
//
// 已知限制：分不出「真值 2ψ」與「ψ 加上大量 50% 次克隆」（兩者在 ψ 尺度下都落在 x.5）。
struct ScaleCheck
{
    bool computed = false;
    double ploidy = 0;           // 判定用的 ψ（4 位小數捨入後）
    double sMax = 0;
    double aMax = 0;
    int segments = 0;
    std::string decision;        // double／keep／unclear
    bool hasCandidate = false;
    FittedPurity candidate;
    bool applied = false;
};

// fits：完整格點（CP-P5 同序）；f0：LOH-TDF 估計（valid 為 false 時不找候選、不改值）。
ScaleCheck checkScale(const std::vector<PurpleCopyNumber> &copyNumbers, const FittedPurity &selected,
        const std::vector<FittedPurity> &fits, bool f0Valid, double f0);

void writeScaleCheckTsv(const std::string &outputDir, const std::string &sampleId, const ScaleCheck &check,
        const FittedPurity &selected);

}

#endif
