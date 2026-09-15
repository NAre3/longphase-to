#ifndef AMBER_SHAREDSCANSINK_H
#define AMBER_SHAREDSCANSINK_H

#include <cstddef>
#include <string>
#include <vector>

#include "BamEvidenceReader.h"

#include "htslib/sam.h"

namespace amber {

// 共用掃描層（EXP-I02）的 AMBER 消費者。
//
// 設計依據 RUN-I002/design.md 的 D1/D2/D9：
//   D1 AMBER 的 RegionTask 只覆蓋 86.1% 基因體，不能反過來當共用層 ⇒ read 由外部送進來。
//   D2 獨立執行下一條 read 會被**所有與它重疊的 region** 各處理一次，
//      因此本 sink 也必須分派給所有重疊的 task，否則 per-locus 計數會少算。
//   D9 一個 sink 只服務**一條染色體**。跨染色體共用單一游標會破壞單調性。
//
// 執行緒模型：每條染色體一個 sink、一個執行緒（`PhasingProcess.cpp` 的 omp 迴圈）。
// RegionTask 對位點是一個分割且不跨染色體，故不同 sink 寫入的 PositionEvidence 互斥，
// 不需要任何鎖。見 design.md「為什麼染色體平行對 AMBER 的 per-locus 計數仍然安全」。
class ContigSink
{
public:
    // tasks 必須是**同一條染色體**、依 start 遞增的連續切片。
    ContigSink(RegionTask *begin, RegionTask *end, int minMappingQuality, int minBaseQuality);

    // record 由共用走訪逐筆送入，必須依座標遞增（BAM 為 coordinate-sorted）。
    void consume(const bam1_t *record);

    const BamScanStats &stats() const { return mStats; }

    // 此 contig 上最大的 task.end。共用走訪的右界取 max(lastSNPpos, 這個值) 即足夠（D9）。
    int maxTaskEnd() const;

private:
    RegionTask *mBegin;
    RegionTask *mEnd;
    std::size_t mFirstCandidate;
    int mMinMappingQuality;
    int mMinBaseQuality;
    BamScanStats mStats;
};

// tasks 依染色體分組的索引。populateTaskQueue 產出的 tasks 天然是「同染色體連續」的
// （它逐一走訪 chrPositions 的分組），本函式只是把每組的 [begin, end) 記下來，
// 並**驗證**該前提成立——若未來 populateTaskQueue 改成交錯輸出，這裡會直接丟例外，
// 不會靜默產生錯誤的計數。
std::vector<std::pair<std::string, std::pair<std::size_t, std::size_t>>>
indexTasksByChromosome(const std::vector<RegionTask> &tasks);

}

#endif
