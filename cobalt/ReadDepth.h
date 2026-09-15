#ifndef COBALT_READDEPTH_H
#define COBALT_READDEPTH_H

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include "../common/ChrBaseRegion.h"

#include "htslib/sam.h"

namespace cobalt {

// DepthReading.java
struct DepthReading
{
    std::string chromosome;      // ChromosomeData.Name（BAM header 原字串，帶 chr 前綴）
    int startPosition = 0;       // windowIndex * WINDOW_SIZE + 1
    double readDepth = 0.0;
    double readGcContent = 0.0;
};

struct ChromosomeSpec
{
    std::string name;
    int length = 0;
};

// ReadDepthAccumulator.java。計數為整數 atomic 加法，可交換 → 與執行緒數無關。
class ReadDepthAccumulator
{
public:
    explicit ReadDepthAccumulator(int windowSize) : mWindowSize(windowSize) {}

    // numWindows = chromosomeLength / windowSize（整數除法，G10）
    void addChromosome(const std::string &chromosome, int chromosomeLength);

    // genomeStart 為 1-based，readStartIndex 為 0-based。thread safe。
    void addReadAlignmentToCounts(const std::string &chromosome, int genomeStart, int alignmentLength,
                                  const std::string &readBases, int readStartIndex);

    std::vector<DepthReading> getChromosomeReadDepths(const std::string &chromosome) const;

private:
    struct ChromosomeWindowCounts
    {
        int numWindows = 0;
        std::unique_ptr<std::atomic<int>[]> baseCounts;
        std::unique_ptr<std::atomic<int>[]> gcCounts;
    };

    const ChromosomeWindowCounts *find(const std::string &chromosome) const;
    ChromosomeWindowCounts *find(const std::string &chromosome);

    int getWindowIndex(int position) const { return (position - 1) / mWindowSize; }
    int getGenomePosition(int windowIndex) const { return windowIndex * mWindowSize + 1; }

    int mWindowSize;
    std::vector<std::string> mNames;
    std::vector<std::unique_ptr<ChromosomeWindowCounts>> mCounts;
};

// BamReadCounter：對每個 partition 平行掃描 BAM，累加到 accumulator。
// 回傳依 chromosomes 順序、window 索引昇冪的 DepthReading（CP-C6 的輸出順序）。
std::vector<DepthReading> calculateReadDepths(const std::string &bamPath,
                                              const std::vector<ChromosomeSpec> &chromosomes,
                                              const std::vector<lp::ChrBaseRegion> &partitions,
                                              int minMappingQuality,
                                              bool includeDuplicates,
                                              int threads);

// ---- 共用掃描層（EXP-I03）----
//
// 整合版不自己開 BAM、不自己建 iterator、不做 partition 路由，read 由 LongPhase-TO 的
// 線性走訪送進來。DepthSink 每條染色體一個，但**共用同一個 accumulator**：
// ReadDepthAccumulator 內部依染色體分槽、各槽自己的 atomic 陣列，因此不同執行緒
// 寫不同染色體的槽互不干擾，不需要鎖。
//
// 為什麼不需要 partition 路由（design.md E1）：partitionGenome 把每條染色體切成
// 連續、不重疊、覆蓋 [1, length] 的區段，而累加是 atomic<int>::fetch_add（可交換）。
// 因此「逐 partition 裁切後加總」與「裁到 [1, length] 一次」結果相同。
class DepthSink
{
public:
    DepthSink(ReadDepthAccumulator &accumulator, std::string chromosome, int contigLength,
              int minMappingQuality, bool includeDuplicates);

    void consume(const bam1_t *record);

private:
    ReadDepthAccumulator &mAccumulator;
    std::string mChromosome;
    int mContigLength;
    int mMinMappingQuality;
    bool mIncludeDuplicates;
};

}

#endif
