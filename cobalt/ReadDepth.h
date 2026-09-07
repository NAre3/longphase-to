#ifndef COBALT_READDEPTH_H
#define COBALT_READDEPTH_H

#include <atomic>
#include <memory>
#include <string>
#include <vector>

#include "../common/ChrBaseRegion.h"

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

}

#endif
