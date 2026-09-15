#ifndef AMBER_BAMEVIDENCEREADER_H
#define AMBER_BAMEVIDENCEREADER_H

#include <cstdint>
#include <string>
#include <vector>

#include "PositionEvidence.h"
#include "../common/SamRecordView.h"

#include "htslib/sam.h"

namespace amber {

// 對應 amber-v4.3 的 RegionTask（RegionTask.java:13-77）。
// Region 的 start/end 為 1-based 含端點，等於該 region 第一個與最後一個位點的座標。
struct RegionTask
{
    std::string chromosome;
    int start = 0;
    int end = 0;
    std::vector<PositionEvidence *> positions;
    std::size_t currentIndex = 0; // 對應 mCurrentIndex，單調前進
};

// 對應 AmberConstants.BAM_MIN_GAP_START（AmberConstants.java:38）
constexpr int BAM_MIN_GAP_START = 4000;

// 對應 AmberConstants 的品質預設值（AmberConstants.java:7-8）
constexpr int DEFAULT_MIN_BASE_QUALITY = 13;
constexpr int DEFAULT_MIN_MAPPING_QUALITY = 50;

// 對應 BamEvidenceReader.populateTaskQueue（BamEvidenceReader.java:67-120）。
// chrPositions 必須已依染色體分組、組內依 position 遞增（對應 TumorAnalysis.java:92-95
// 的 Collections.sort，該排序在切分之前執行，切分規則依序走訪這個 List）。
// 切分條件為嚴格小於：end + minGap == pos 併入同一個 region。
std::vector<RegionTask> populateTaskQueue(
        const std::vector<std::pair<std::string, std::vector<PositionEvidence *>>> &chrPositions,
        int minGap);

struct BamScanStats
{
    std::uint64_t recordsConsumed = 0;   // 通過 BamSlicerFilter 而進入 processRecord 的 read 數
    std::uint64_t nonAcgtnBases = 0;     // 取到非 ACGTN 鹼基碼的次數（Java 端會在此丟例外）
};

// 對應 BamEvidenceReader.processBam + BamReaderThread + PositionEvidenceChecker。
//
// threads > 1 時以 region 為單位分工。這在結構上是安全的：region 對位點是一個分割
// （populateTaskQueue 把每個位點恰好放進一個 region），因此不同執行緒寫入的
// PositionEvidence 互不重疊；每個 region 內部的 read 走訪順序仍由 htslib 的
// iterator 決定，與執行緒數無關。**但這是推論，須以實測確認**——
// Java 端的對應結論由 EXP-002 實測建立（threads 1 vs 8，12 個 checkpoint 逐位元組相同）。
BamScanStats processBam(
        const std::string &bamFile, std::vector<RegionTask> &tasks,
        int minMappingQuality, int minBaseQuality, int threads);

// ---- 共用掃描層（EXP-I02）用的兩個入口 ----
//
// 整合版不自己開 BAM、不自己建 iterator，read 由 LongPhase-TO 的線性走訪送進來，
// 因此需要把 processBam 內部的兩段邏輯轉出來給 ContigSink 用。
// 兩者都只是轉呼叫既有實作，**沒有任何行為分支**——這是「integrated 與 amber_port
// 走同一份程式碼」這句話的依據。

// 對應 BamSlicerFilter.passesFilters：0x4 / 0x100 / 0x800 / 0x400 任一命中即排除。
// **不含 MAPQ**——AMBER 的 MAPQ 是在 addEvidence 內計數而非丟棄（RUN-I001 行為對照表）。
bool passesSlicerFilters(const bam1_t *record);

// 對應 RegionTask.processRecord。task.currentIndex 單調前進。
void processRecordForRegion(RegionTask &task, const lp::SamRecordView &read,
        int minMappingQuality, int minBaseQuality, BamScanStats &stats);

}

#endif
