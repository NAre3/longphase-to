#ifndef AMBER_TUMORFILTERS_H
#define AMBER_TUMORFILTERS_H

#include <string>
#include <vector>

#include "PositionEvidence.h"

namespace amber {

// 對應 AmberConstants 的 tumor-only 預設值（AmberConstants.java:9,12,27）。
// TumorMinDepth 在 tumor-only（ReferenceIds 為空）下取 DEFAULT_TUMOR_ONLY_MIN_DEPTH 25，
// 而非 tumor/normal 的 DEFAULT_TUMOR_MIN_DEPTH 8（AmberConfig.java:135-146）。
constexpr int DEFAULT_TUMOR_ONLY_MIN_DEPTH = 25;
constexpr int DEFAULT_TUMOR_ONLY_MIN_SUPPORT = 2;
constexpr double QUAL_FILTERED_THRESHOLD = 0.15;

// 對應 AmberUtils.aboveQualFilter（AmberUtils.java:48-56）。
// 分母是 ReadDepth；三個 filtered 計數器相加後取比例，**嚴格小於**門檻才通過。
bool aboveQualFilter(const PositionEvidence &pe);

// 對應 HumanChromosome.chromosomeRank（HumanChromosome.java:162-190）：
// 去 chr 前綴後 1-22 取其數值、X 為 23、Y 為 24、MT/M 為 25，其餘為 26。
int chromosomeRank(const std::string &chromosome);

// 對應 GenomePosition.compare（GenomePosition.java:31-39）：
// 先以 ContigComparator 比 rank（ContigComparator.java:10-19），rank 相同再比 position。
// 注意 rank 是**數值序**，不是字串序——chr2 排在 chr10 之前。
bool genomePositionLess(const PositionEvidence &a, const PositionEvidence &b);

}

#endif
