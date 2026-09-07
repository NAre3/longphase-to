#include "Segmentation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>

#include "CobaltConstants.h"
#include "../common/HumanChromosome.h"

namespace cobalt {

namespace {

// CobaltRatioSegmenter.buildSegmentationData（第 80–103 行）
//   rawValues[i] = v            （在 floor 判斷之前）
//   v < 0.001 -> -9.965784      （硬編碼 floor，嚴格小於）
//   否則      -> (float) FastMath.log(2, v) = (float)(log(v)/log(2))
// **(float) 轉型屬演算法邏輯必須保留**（G4）——它是 DECISION D2 的成因。
// FastMath.log(base, x) = log(x)/log(base)，**base 在前**，寫反會得到 log(2)/log(v)。
inline double valueForSegmentation(double v)
{
    if(v < 0.001){ return -9.965784; }
    return static_cast<double>(static_cast<float>(std::log(v) / std::log(2.0)));
}

// Doubles.mean：單純循序累加除以長度（位元組碼確認）
double doublesMean(const double *values, std::size_t count)
{
    if(count == 0){ return std::nan(""); }
    double sum = 0.0;
    for(std::size_t i = 0; i < count; ++i){ sum += values[i]; }
    return sum / static_cast<double>(count);
}

// CobaltRatioFile.FORMAT = DecimalFormat("#.####")：HALF_EVEN + 去尾隨零
std::string format4(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4f", v);
    std::string s(buf);
    if(s.find('.') != std::string::npos)
    {
        while(!s.empty() && s.back() == '0'){ s.pop_back(); }
        if(!s.empty() && s.back() == '.'){ s.pop_back(); }
    }
    if(s == "-0"){ s = "0"; }
    return s;
}

}

SegmentationResult segmentRatios(const std::vector<CobaltRatio> &ratios, double gamma)
{
    SegmentationResult result;
    result.gamma = gamma;
    result.isWindowed = true;                       // CobaltRatioSegmenter.isWindowed()（G3）

    // ---- PerArmSegmenter 建構：value(r) = tumorGCRatio，過濾 >= 0.0（非嚴格）----
    std::map<std::string, std::size_t> armIndex;
    for(const CobaltRatio &r : ratios)
    {
        const double v = r.tumorGCRatio;
        if(!(v >= 0.0)){ continue; }

        const std::string shortName = lp::stripChrPrefix(r.chromosome);
        const char arm = r.position < lp::centromereOf(r.chromosome) ? 'P' : 'Q';
        const std::string armId = shortName + "_" + arm;

        auto found = armIndex.find(armId);
        if(found == armIndex.end())
        {
            armIndex.emplace(armId, result.arms.size());
            ArmData a;
            a.armId = armId;
            a.chromosomeShort = shortName;
            a.chromosomeRank = lp::chromosomeOrdinal(shortName);
            a.arm = arm;
            result.arms.push_back(std::move(a));
            found = armIndex.find(armId);
        }
        ArmData &a = result.arms[found->second];
        a.rawValues.push_back(v);                                   // 在 floor 之前
        a.valuesForSegmentation.push_back(valueForSegmentation(v));
        a.positions.push_back(r.position);
    }

    // ChrArm.compareTo：先比染色體 enum 序，再比 arm（P < Q）
    std::sort(result.arms.begin(), result.arms.end(), [](const ArmData &a, const ArmData &b){
        if(a.chromosomeRank != b.chromosomeRank){ return a.chromosomeRank < b.chromosomeRank; }
        return a.arm < b.arm;
    });

    int totalCount = 0;
    for(const ArmData &a : result.arms){ totalCount += static_cast<int>(a.valuesForSegmentation.size()); }
    result.totalCount = totalCount;
    result.penaltyMode = totalCount < UNIFORM_PENALTY_THRESHOLD ? "uniform" : "per-arm-gamma";

    if(result.penaltyMode == "uniform")
    {
        // 比照 amber/Segmentation.cpp 的處置：allRatios 的組裝順序取自 mDataByArm.keySet()
        // （identity hash），跨語言不可重現。寧可明確失敗，也不要用未驗證的順序默默算下去。
        throw std::runtime_error(
            "uniform penalty branch is not implemented: allRatios assembly order comes from "
            "mDataByArm.keySet() and is not reproducible across languages. totalCount="
            + std::to_string(totalCount));
    }

    // ---- 逐 arm：penalty -> 分段 -> segment 座標 ----
    for(ArmData &a : result.arms)
    {
        if(a.valuesForSegmentation.empty()){ continue; }

        const double penalty = lp::gammaSegmentPenalty(a.valuesForSegmentation, gamma, true, &a.trace);
        a.fit = lp::segment(a.valuesForSegmentation, penalty, &a.pcfMeans);

        std::size_t idx = 0;
        for(std::size_t i = 0; i < a.fit.lengths.size(); ++i)
        {
            const int count = a.fit.lengths[i];
            PcfSegmentOut seg;
            seg.chromosome = a.chromosomeShort;
            seg.start = a.positions[idx];
            // ChromosomeArmSegments.MeanRatio 取 **rawValues** 該段的平均（非 pcfMeans）
            seg.meanRatio = doublesMean(a.rawValues.data() + idx, static_cast<std::size_t>(count));
            idx += static_cast<std::size_t>(count);
            // WindowSegments.segmentEnd = endRatio.position() + WINDOW_SIZE - 1（G3）
            seg.end = a.positions[idx - 1] + WINDOW_SIZE - 1;
            a.segments.push_back(std::move(seg));
        }
    }

    return result;
}

void writeSegmentsFile(const std::string &path, const SegmentationResult &result)
{
    std::ofstream out(path);
    if(!out){ throw std::runtime_error("cannot open for write: " + path); }
    // SegmentsFile.java:25 的硬編碼表頭，**四欄，無 n.probes**（findings F1）
    out << "Chromosome\tStart\tEnd\tMeanRatio\n";
    for(const ArmData &a : result.arms)
    {
        for(const PcfSegmentOut &s : a.segments)
        {
            out << "chr" << s.chromosome << '\t' << s.start << '\t' << s.end << '\t'
                << format4(s.meanRatio) << '\n';
        }
    }
}

}
