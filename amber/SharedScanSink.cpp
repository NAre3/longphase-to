#include "SharedScanSink.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

#include "../common/SamRecordView.h"

namespace amber {

using namespace lp;

ContigSink::ContigSink(RegionTask *begin, RegionTask *end, int minMappingQuality, int minBaseQuality)
    : mBegin(begin), mEnd(end), mFirstCandidate(0),
      mMinMappingQuality(minMappingQuality), mMinBaseQuality(minBaseQuality)
{
}

int ContigSink::maxTaskEnd() const
{
    int maxEnd = 0;
    for(RegionTask *task = mBegin; task != mEnd; ++task){
        maxEnd = std::max(maxEnd, task->end);
    }
    return maxEnd;
}

void ContigSink::consume(const bam1_t *record)
{
    // 共用層零過濾（RUN-I001 的結論），AMBER 的 slicer filter 在此消費者內執行。
    if(!passesSlicerFilters(record)){
        return;
    }

    // recordsConsumed：**刻意不重現**獨立執行下的 per-region 造訪次數語義（design.md D3，
    // 使用者 2026-09-15 裁示）。此處計的是「通過 slicer filter 的 read 數」，每條 read 一次。
    // 該欄位的唯一消費者是一行 stderr 日誌，不進任何 checkpoint 或 stage 輸出。
    ++mStats.recordsConsumed;

    const SamRecordView read(record);
    const int alignmentStart = read.alignmentStart();
    const int alignmentEnd = read.alignmentEnd();

    const std::size_t taskCount = static_cast<std::size_t>(mEnd - mBegin);

    // 單調前進：整段落在這條 read 左側的 task 之後都不會再被用到。
    // 成立的前提是 read 依座標遞增抵達，且 task 依 start 遞增排列。
    while(mFirstCandidate < taskCount && mBegin[mFirstCandidate].end < alignmentStart){
        ++mFirstCandidate;
    }

    // 分派給**所有**重疊的 task（D2）。task 依 start 遞增，故第一個 start > alignmentEnd
    // 之後的都不重疊，可提前收手。
    for(std::size_t index = mFirstCandidate; index < taskCount; ++index){
        RegionTask &task = mBegin[index];

        if(task.start > alignmentEnd){
            break;
        }

        // mFirstCandidate 之後仍可能有 end < alignmentStart 的 task
        // （task 之間長度不一，區間不是巢狀的），故逐一再檢查一次。
        if(task.end < alignmentStart){
            continue;
        }

        processRecordForRegion(task, read, mMinMappingQuality, mMinBaseQuality, mStats);
    }
}

std::vector<std::pair<std::string, std::pair<std::size_t, std::size_t>>>
indexTasksByChromosome(const std::vector<RegionTask> &tasks)
{
    std::vector<std::pair<std::string, std::pair<std::size_t, std::size_t>>> index;
    std::unordered_set<std::string> seen;

    std::size_t i = 0;
    while(i < tasks.size()){
        const std::string &chromosome = tasks[i].chromosome;

        if(!seen.insert(chromosome).second){
            // populateTaskQueue 逐組輸出，同一條染色體不該出現第二個區段。
            throw std::runtime_error(
                    "amber tasks are not grouped by chromosome: " + chromosome);
        }

        const std::size_t begin = i;
        while(i < tasks.size() && tasks[i].chromosome == chromosome){
            if(i > begin && tasks[i].start < tasks[i - 1].start){
                throw std::runtime_error(
                        "amber tasks are not sorted by start on " + chromosome);
            }
            ++i;
        }

        index.emplace_back(chromosome, std::make_pair(begin, i));
    }

    return index;
}

}
