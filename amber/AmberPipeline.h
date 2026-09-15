#ifndef AMBER_AMBERPIPELINE_H
#define AMBER_AMBERPIPELINE_H

#include <string>
#include <vector>

#include "BamEvidenceReader.h"
#include "PositionEvidence.h"

namespace amber {

// AMBER tumor-only 流程的參數。欄位與 amber_port 的 CLI 一一對應
// （AmberApplication.cpp 的 usage()），預設值與該處相同。
struct PipelineConfig
{
    std::string lociPath;
    std::string bedPath;
    std::string tumorBam;
    std::string debugOnlyChr;
    std::string outputDir;
    std::string sampleId;
    int threads = 1;
    bool writeTumorData = false;
    bool writeVersion = false;
    int minBaseQuality = DEFAULT_MIN_BASE_QUALITY;
    int minMappingQuality = DEFAULT_MIN_MAPPING_QUALITY;
};

// prescan 的產出。tasks 內的 PositionEvidence* 指向 evidence 的元素，
// 因此 **evidence 在 postscan 之前不得再被 push_back 或複製**
// （move 是安全的：std::vector 搬移時堆積緩衝區不變）。
struct PrescanResult
{
    std::vector<PositionEvidence> evidence;
    std::vector<RegionTask> tasks;
    // tumorBam 未給定時 amber_port 只跑到 CP-A2 就結束；此旗標保留該行為。
    bool hasBamStage = false;
};

// ---- 流程的三段 ----
//
// 拆點由「BAM 掃描發生在哪裡」決定，這是 EXP-I02 的整合要求：
//   prescan   CP-A1 → CP-A2 → CP-A2b（site 載入、blacklist、region 切分）—— 掃描之前
//   （掃描）  amber_port 走 processBam；整合版走共用掃描層的 ContigSink
//   postscan  CP-A3 → CP-A9 與三個 stage 輸出 —— 掃描之後
//
// amber_port 的 main 是這三段的薄包裝，**與整合版呼叫同一組函式**。
// 這是「凍結候選的行為未因整合而改變」這句話的依據。
PrescanResult prescan(const PipelineConfig &cfg);

void postscan(const PipelineConfig &cfg, const std::vector<PositionEvidence> &evidence,
        const BamScanStats &stats);

}

#endif
