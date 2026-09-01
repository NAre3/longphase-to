#ifndef AMBER_BAMEVIDENCEREADER_H
#define AMBER_BAMEVIDENCEREADER_H

#include <cstdint>
#include <string>
#include <vector>

#include "PositionEvidence.h"

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
// 單執行緒；RUN-002 已於全基因體確認 threads 1 與 8 的 12 個 checkpoint 逐位元組相同。
BamScanStats processBam(
        const std::string &bamFile, std::vector<RegionTask> &tasks,
        int minMappingQuality, int minBaseQuality);

}

#endif
