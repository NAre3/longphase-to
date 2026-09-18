#ifndef COBALT_COBALTPIPELINE_H
#define COBALT_COBALTPIPELINE_H

#include <string>
#include <vector>

#include "../common/ChrBaseRegion.h"
#include "CobaltOutput.h"
#include "CobaltWindow.h"
#include "Segmentation.h"
#include "GcProfile.h"
#include "ReadDepth.h"
#include "Regions.h"
#include "WindowStatuses.h"

namespace cobalt {

// COBALT tumor-only 流程的參數。欄位與 cobalt_port 的 CLI 一一對應
// （CobaltApplication.cpp 的 usage），預設值與該處相同。
struct PipelineConfig
{
    std::string bamPath;
    std::string gcProfile;
    std::string diploidBed;
    std::string excludedPath;
    std::string outputDir = ".";
    std::string sampleId = "tumor";
    int threads = 1;
    int minMappingQuality = 10;      // cobalt_port 的 -min_quality 預設
    double pcfGamma = 100.0;         // -pcf_gamma 預設；AMBER/COBALT 兩邊皆硬編碼 100
    bool includeDuplicates = false;
    // 三個 stage 輸出是否落檔。**預設 true**：既有呼叫端（cobalt_port、
    // longphase-to 的 EXP-I03 整合）行為完全不變。
    // 設 false 只在下游改以記憶體接手結果時使用（PURPLE 整合）。
    bool writeStageOutputs = true;
};

// prescan 的產出（CP-C1 → CP-C5）。postscan 需要其中四項：
// chromosomes（決定 CP-C6 的輸出順序）、gcData、statuses，以及 partitions
// ——partitions 只有 cobalt_port 的獨立掃描會用到，整合版不用（design.md E1）。
struct PrescanResult
{
    std::vector<ChromosomeSpec> chromosomes;
    std::vector<lp::ChrBaseRegion> partitions;
    GcProfileData gcData;
    std::vector<ExcludedRegion> excluded;
    DiploidRegions diploid;
    WindowStatuses statuses;
};

// ---- 流程的三段 ----
//
// 拆點由「BAM 掃描發生在哪裡」決定（EXP-I03）：
//   prescan   CP-C1 → CP-C5（染色體、partition、GC profile、excluded/diploid、window status）
//   （掃描）  cobalt_port 走 calculateReadDepths；整合版走共用掃描層的 DepthSink
//   postscan  CP-C6 → CP-C14 與三個 stage 輸出
//
// cobalt_port 的 main 是這三段的薄包裝，**與整合版呼叫同一組函式**。
// 這是「凍結候選（b71e68e）的行為未因整合而改變」這句話的依據。
PrescanResult prescan(const PipelineConfig &cfg);

// postscan 的產出。這兩項是 PURPLE 需要的兩個 stage 輸出的記憶體來源：
//   ratios        -> <sample>.cobalt.ratio.tsv.gz
//   segmentation  -> <sample>.cobalt.ratio.pcf
// （gc.median 只有 COBALT 自己的診斷用途，PURPLE 不讀，故不回傳。）
struct PostscanResult
{
    std::vector<CobaltRatio> ratios;
    SegmentationResult segmentation;
};

PostscanResult postscan(const PipelineConfig &cfg, const PrescanResult &pre,
        const std::vector<DepthReading> &depths);

}

#endif
