#include "Segmentation.h"
#include "../common/Segmentation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace amber {

// 通用分段核心已抽到 common/（namespace lp），此處只留 AMBER 形狀的 wrapper。
using namespace lp;


SegmentationResult segmentBafs(const std::vector<AmberBAF> &bafs)
{
    SegmentationResult result;

    // PerArmSegmenter 建構（PerArmSegmenter.java:39-52）：
    // 走訪各染色體的 BAF，value(baf) = tumorModifiedBAF，>= 0.0 者依 (染色體, arm) 分組。
    // tumorModifiedBAF = 0.5 + |baf - 0.5| >= 0.5，故實務上全部入選。
    std::map<std::string, std::size_t> armIndex;

    for(const AmberBAF &baf : bafs){
        const double value = baf.tumorModifiedBAF();
        if(!(value >= 0.0)){
            continue;
        }

        const std::string shortName = stripChrPrefix(baf.chromosome);
        const char arm = baf.position < centromereOf(baf.chromosome) ? 'P' : 'Q';
        const std::string armId = shortName + "_" + arm;

        auto found = armIndex.find(armId);
        if(found == armIndex.end()){
            armIndex.emplace(armId, result.arms.size());
            ArmSegments segments;
            segments.armId = armId;
            segments.chromosomeShort = shortName;
            segments.chromosomeRank = chromosomeOrdinal(shortName);
            segments.arm = arm;
            result.arms.push_back(std::move(segments));
            found = armIndex.find(armId);
        }

        result.arms[found->second].values.push_back(value);
    }

    // ChrArm.compareTo：先比染色體 enum 序，再比 arm（P < Q）
    std::sort(result.arms.begin(), result.arms.end(), [](const ArmSegments &a, const ArmSegments &b){
        if(a.chromosomeRank != b.chromosomeRank){
            return a.chromosomeRank < b.chromosomeRank;
        }
        return a.arm < b.arm;
    });

    int totalCount = 0;
    for(const ArmSegments &arm : result.arms){
        totalCount += static_cast<int>(arm.values.size());
    }
    result.totalCount = totalCount;
    result.penaltyMode = totalCount < UNIFORM_PENALTY_THRESHOLD ? "uniform" : "per-arm-gamma";

    if(result.penaltyMode == "uniform"){
        // PerArmSegmenter.java:55-72 的 uniform 分支。其 allRatios 的組裝順序取自
        // HashMap.keySet()，跨語言不保證。本研究的資料不走此分支（spec U5 已關閉），
        // 故不實作——寧可明確失敗，也不要用一個未經驗證的順序默默算出結果。
        throw std::runtime_error(
                "uniform penalty branch is not implemented: its allRatios assembly order comes from "
                "HashMap.keySet() and is not reproducible across languages (see spec U5). "
                "totalCount=" + std::to_string(totalCount));
    }

    // 需要每個 arm 內的 BAF（依原順序）以取得 segment 的起訖座標
    std::map<std::string, std::vector<const AmberBAF *>> bafsByArm;
    for(const AmberBAF &baf : bafs){
        if(!(baf.tumorModifiedBAF() >= 0.0)){
            continue;
        }
        const std::string shortName = stripChrPrefix(baf.chromosome);
        const char arm = baf.position < centromereOf(baf.chromosome) ? 'P' : 'Q';
        bafsByArm[shortName + "_" + arm].push_back(&baf);
    }

    for(ArmSegments &arm : result.arms){
        if(arm.values.empty()){
            continue;
        }

        // Segmenter 以 GammaPenaltyCalculator(gamma, true) 逐 arm 計算 penalty
        const double penalty = gammaSegmentPenalty(arm.values, BAF_SEGMENTATION_GAMMA, true);
        const Fit fit = segment(arm.values, penalty);

        // ChromosomeArmSegments 建構（ChromosomeArmSegments.java:14-30）：
        // MeanRatio 取 rawValues 該段的平均（Doubles.mean），**不是** PiecewiseConstantFit
        // 裡那個已四捨五入到三位小數的 means。
        const std::vector<const AmberBAF *> &armBafs = bafsByArm[arm.armId];
        std::size_t ratiosIndex = 0;

        for(std::size_t i = 0; i < fit.lengths.size(); ++i){
            const int count = fit.lengths[i];
            const double meanRatio = mean(arm.values.data() + ratiosIndex, static_cast<std::size_t>(count));

            const AmberBAF *startRatio = armBafs[ratiosIndex];
            ratiosIndex += static_cast<std::size_t>(count);
            const AmberBAF *endRatio = armBafs[ratiosIndex - 1];

            PcfSegment segment;
            segment.chromosome = arm.chromosomeShort;
            segment.start = startRatio->position;
            segment.end = endRatio->position; // AbsoluteSegments.segmentEnd = endRatio.position()
            segment.meanRatio = meanRatio;
            arm.segments.push_back(std::move(segment));
        }
    }

    return result;
}

namespace {

// 對應 CobaltRatioFile.FORMAT = DecimalFormat("#.####")：
// 最多四位小數、去掉尾端的 0 與孤立的小數點。
std::string format4Significant(double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);
    std::string s = buffer;

    if(s.find('.') != std::string::npos){
        while(!s.empty() && s.back() == '0'){
            s.pop_back();
        }
        if(!s.empty() && s.back() == '.'){
            s.pop_back();
        }
    }

    return s;
}

}

void writeSegmentsFile(const std::string &path, const SegmentationResult &result)
{
    std::ofstream out(path);
    if(!out){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    out << "Chromosome\tStart\tEnd\tMeanRatio\n";

    for(const ArmSegments &arm : result.arms){
        for(const PcfSegment &segment : arm.segments){
            // SegmentsFile 以 genomeVersion.versionedChromosome 加上 chr 前綴
            out << "chr" << arm.chromosomeShort << '\t'
                << segment.start << '\t'
                << segment.end << '\t'
                << format4Significant(segment.meanRatio) << '\n';
        }
    }
}


void writeSegmentationCheckpoints(const SegmentationResult &result)
{
    if(!CpDump::enabled()){
        return;
    }

    // CP-A8：分段輸入與 penalty（對應 patch_segmenter.py 插在 PerArmSegmenter 建構子尾端）
    std::vector<CpDump::Row> a8;
    a8.push_back({"", 0, "totalCount\t" + std::to_string(result.totalCount)});
    a8.push_back({"", 0, "uniformPenaltyThreshold\t" + std::to_string(result.uniformPenaltyThreshold)});
    a8.push_back({"", 0, "penaltyMode\t" + result.penaltyMode});
    a8.push_back({"", 0, "gamma\t" + CpDump::num(result.gamma)});

    for(const ArmSegments &arm : result.arms){
        a8.push_back({"", 0, "armCount." + arm.armId + "\t" + std::to_string(arm.values.size())});
        for(std::size_t i = 0; i < arm.values.size(); ++i){
            // BAFSegmenter 的 valuesForSegmentation 與 rawValues 相同（BAFSegmenter.java:47-58）
            a8.push_back({"", 0, "armValue." + arm.armId + "." + std::to_string(i) + "\t"
                    + CpDump::num(arm.values[i]) + "\t" + CpDump::num(arm.values[i])});
        }
    }
    CpDump::write("CP-A8", "field\tvalue\trawValue", a8);

    // CP-A9：分段結果
    std::vector<CpDump::Row> a9;
    for(const ArmSegments &arm : result.arms){
        int idx = 0;
        for(const PcfSegment &segment : arm.segments){
            a9.push_back({"", 0, arm.armId + "\t" + std::to_string(idx++) + "\t"
                    + segment.chromosome + "\t" + std::to_string(segment.start) + "\t"
                    + std::to_string(segment.end) + "\t" + CpDump::num(segment.meanRatio)});
        }
    }
    CpDump::write("CP-A9", "arm\tidx\tchromosome\tstart\tend\tmeanRatio", a9);
}

}
