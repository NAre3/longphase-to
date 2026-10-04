#include "PurpleInputAdapter.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "PurpleInput.h"

namespace purple {

namespace {

// amber/AmberOutput.cpp:19 format4 的複本。
std::string amberFormat4(double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);
    return buffer;
}

// cobalt/CobaltOutput.cpp:17 decimalFormat4 的複本。
std::string cobaltDecimalFormat4(double v)
{
    if(std::isnan(v)){ return "NaN"; }
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

double roundTripAmberBaf(double value)
{
    // AmberBAFFile 對非有限值寫字面 "0"（AmberOutput.cpp:30 formatOrZero），
    // 讀回即 0.0。這條路徑必須一併重現，否則 NaN/Inf 會以原值進入 PURPLE。
    if(!std::isfinite(value)){ return 0.0; }
    return std::stod(amberFormat4(value));
}

double roundTripCobaltRatio(double value)
{
    const std::string text = cobaltDecimalFormat4(value);
    // "NaN" 經 std::stod 得到 NaN，與 PurpleInput 讀檔時的結果一致。
    return std::stod(text);
}

InputData buildInputsFromMemory(
        const std::string &sampleId,
        const std::vector<amber::AmberBAF> &bafs,
        const amber::SegmentationResult &amberSegmentation,
        const std::vector<cobalt::CobaltRatio> &ratios,
        const cobalt::SegmentationResult &cobaltSegmentation,
        double contamination)
{
    InputData data;
    data.sampleId = sampleId;
    data.contamination = contamination;

    // ---- BAF ----
    // 欄位對應 PurpleInput.cpp:88-95 讀回來的那六欄；深度是整數，不需捨入。
    data.bafs.reserve(bafs.size());
    for(const amber::AmberBAF &in : bafs){
        AmberBaf out;
        out.chromosome = in.chromosome;
        out.position = in.position;
        out.tumorBaf = roundTripAmberBaf(in.tumorBAF);
        out.tumorDepth = in.tumorDepth;
        out.normalBaf = roundTripAmberBaf(in.normalBAF);
        out.normalDepth = in.normalDepth;
        data.bafs.push_back(std::move(out));
    }

    // gender 與平均深度由 BAF 推導，**必須在 BAF 捨入之後算**，
    // 否則與讀檔路徑的輸入不同（PurpleInput.cpp:218-220 是在讀檔後才算）。
    data.averageTumorDepth = averageTumorDepth(data.bafs);
    data.amberGender = determineAmberGender38(data.bafs);
    data.cobaltGender = data.amberGender;

    // ---- PCF ----
    // 兩個寫檔端都以 "chr" + 短名輸出（amber/Segmentation.cpp:162、
    // cobalt/Segmentation.cpp:148），而 PURPLE 讀檔時 end 會 +1
    // （PurpleInput.cpp 新格式分支）。此處重現同樣的換算，
    // 並依寫檔端的 arm/segment 走訪順序產生列序。
    std::vector<PcfInterval> amberIntervals;
    for(const amber::ArmSegments &arm : amberSegmentation.arms){
        for(const amber::PcfSegment &segment : arm.segments){
            amberIntervals.push_back(PcfInterval{
                    "chr" + arm.chromosomeShort, segment.start, segment.end + 1});
        }
    }
    data.amberPcf = buildPcfPositions(amberIntervals, PcfSource::TUMOR_BAF);

    std::vector<PcfInterval> cobaltIntervals;
    for(const cobalt::ArmData &arm : cobaltSegmentation.arms){
        for(const cobalt::PcfSegmentOut &segment : arm.segments){
            cobaltIntervals.push_back(PcfInterval{
                    "chr" + segment.chromosome, segment.start, segment.end + 1});
        }
    }
    data.cobaltTumorPcf = buildPcfPositions(cobaltIntervals, PcfSource::TUMOR_RATIO);

    // ---- ratio ----
    // 除了四位小數捨入，還要重現 PurpleInput.cpp:101-127 讀檔時做的兩件事：
    // referenceReadDepth == -1 時把兩個 reference ratio 視為 1，以及
    // 對 reference ratio 套 genderAdjusted。後者在讀檔路徑是用
    // 讀檔當下的 amberGender，這裡用的是同一個值。
    data.ratios.reserve(ratios.size());
    for(const cobalt::CobaltRatio &in : ratios){
        data.ratios.push_back(makeCobaltRatio(
                in.chromosome,
                in.position,
                roundTripCobaltRatio(in.referenceReadDepth),
                roundTripCobaltRatio(in.referenceGCRatio),
                roundTripCobaltRatio(in.referenceGCDiploidRatio),
                roundTripCobaltRatio(in.referenceGCContent),
                roundTripCobaltRatio(in.tumorReadDepth),
                roundTripCobaltRatio(in.tumorGCRatio),
                roundTripCobaltRatio(in.tumorGcContent),
                data.amberGender));
    }

    return data;
}

}
