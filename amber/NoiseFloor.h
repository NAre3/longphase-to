#ifndef AMBER_NOISEFLOOR_H
#define AMBER_NOISEFLOOR_H

#include <string>
#include <vector>

#include "PositionEvidence.h"

namespace amber {

// 對應 amber-v4.3 的 amber.purity 套件（TumorOnlyPurityAnalysis 與其相依類別）。
// 行為出處逐條見 NoiseFloor.cpp 的註解。

// AmberConstants.java:43-58 的常數
constexpr double LOWER_CDF_BOUND_FOR_CAPTURE = 0.16;
constexpr double UPPER_CDF_BOUND_FOR_CAPTURE = 0.84;
constexpr int MINIMUM_CAPTURED_POINTS = 15;
constexpr double PEAK_SEARCH_START = 0.005;
constexpr double PEAK_SEARCH_END = 0.37;
constexpr double PEAK_SEARCH_STEP_RATIO = 1.05;
constexpr double PEAK_SEARCH_INITIAL_STEP = 0.001;
constexpr double PEAK_SEARCH_OVERSHOOT = 0.0001;
constexpr double HET_VAF_LOWER_BOUND = 0.35;
constexpr double HET_VAF_UPPER_BOUND = 0.65;
constexpr double MIN_CUTOFF = 0.04;
constexpr double HOMOZYGOUS_PROPORTION_LOWER_BOUND_FOR_CONTAMINATION = 0.25;
constexpr double HOMOZYGOUS_PROPORTION_UPPER_BOUND_FOR_CONTAMINATION = 0.75;
constexpr double MUTATION_AUC_LOWER_BOUND_FOR_CONTAMINATION = 0.825;
constexpr double CHR_ARM_LOWER_BOUND_FOR_CONTAMINATION = 0.5;
constexpr double CHR_ARM_AUC_UPPER_BOUND_FOR_COPY_NUMBER_EVENTS = 0.1;

// Doubles.java:18,45 —— greaterOrEqual 是 epsilon 比較，不是 >=
constexpr double DOUBLES_EPSILON = 1e-10;
inline bool doublesGreaterOrEqual(double value, double reference)
{
    return value - reference > -DOUBLES_EPSILON;
}
inline bool doublesIsZero(double value)
{
    double diff = value - 0.0;
    return (diff < 0 ? -diff : diff) < DOUBLES_EPSILON;
}

struct PeakDiagnostic
{
    double vaf;
    double score;
    double homozygousProportion;
    double chrArmAuc;
    double mutationAuc;
    int capturedPoints;
    std::string classification; // "contamination" / "copy-number" / "none" / "gnomad-rejected"
};

struct NoiseFloorResult
{
    double noiseFloor = 0.0;
    double contamination = 0.0;
    std::vector<double> contaminationPeakVafs;

    // 診斷用（不進驗收規則）：Java 端以 -log_debug 印出同樣的量
    std::size_t evidencePoints = 0;
    std::size_t evidencePointsAfterImmuneFilter = 0;
    double baselineHetGnomadFrequency = 0.0;
    std::vector<PeakDiagnostic> maximaDiagnostics;
    std::vector<std::pair<double, double>> gridScores; // 每個 grid level 的 (vaf, score)
};

// 對應 TumorOnlyPurityAnalysis 的建構與 cutoff()。
// evidence 為 CP-A5 的 rawData（四道 filter 與排序之後）。
NoiseFloorResult computeNoiseFloor(const std::vector<PositionEvidence> &evidence);

}

#endif
