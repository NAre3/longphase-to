#include "BamEvidenceReader.h"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <thread>

#include "../common/SamRecordView.h"

#include "htslib/hts.h"
#include "htslib/sam.h"

namespace amber {

using namespace lp;

namespace {

// 對應 BamSlicerFilter.passesFilters，參數為 BamSlicer(0, false, false, false)
// （BamReaderThread.java:109；BamSlicerFilter.java:35-61）。
// minMappingQuality 為 0 時 `getMappingQuality() < 0` 恆為假 → **此層不做 MAPQ 過濾**。
// MAPQ 的門檻在 addEvidence 內，而且是計數而非丟棄。
constexpr uint16_t SLICER_EXCLUDED_FLAGS =
        BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FDUP;

inline bool passesFilters(const bam1_t *record)
{
    return (record->core.flag & SLICER_EXCLUDED_FLAGS) == 0;
}

inline bool isAcgtn(char base)
{
    return base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N';
}

// 對應 PositionEvidenceChecker.getBaseQuality（PositionEvidenceChecker.java:99-113）。
// 位點落在 deletion / ref skip 上時往右掃到第一個有 read position 的參考位置；
// 掃到 alignmentEnd 仍無則回傳 0（0 < MinBaseQuality → 計入 baseQualFiltered）。
int getBaseQuality(int position, const SamRecordView &read)
{
    for(int pos = position; pos <= read.alignmentEnd(); ++pos){
        const int readPosition = read.readPositionAtReferencePosition(pos);

        if(readPosition > 0){
            return read.baseQualityAt(readPosition);
        }
    }

    return 0;
}

// 對應 PositionEvidenceChecker.isIndel（PositionEvidenceChecker.java:82-97）。
// readIndex 參數為 1-based read position，與 Java 的呼叫方式一致（該處傳入 readIndex + 1，
// 註解寫明 "to preserve existing read index behaviour"）。+1/-1 的堆疊刻意不化簡。
bool isIndel(int bafPosition, int readIndex, const SamRecordView &read)
{
    if(read.alignmentEnd() > bafPosition){
        // Delete?
        if(read.readPositionAtReferencePosition(bafPosition + 1) == 0){
            return true;
        }

        // Insert?
        return read.referencePositionAtReadPosition(readIndex + 1) != bafPosition + 1;
    }

    return false;
}

// 對應 PositionEvidenceChecker.addEvidence（PositionEvidenceChecker.java:23-80）。
// 必須逐條保留的非直覺行為見 behaviour-contract.md §2.1：
//   - readDepth 在 filter 判定之後、return 之前無條件遞增（被擋掉的 read 仍計入）
//   - 兩個品質 filter 各自獨立判定，同一條 read 可同時計入兩個計數器
//   - readIndex < 0（位點落在 deletion）時 Ref/Alt 與 indelCount 三者都不動
//   - seqTechFiltered 只在 ULTIMA 下遞增，本研究恆為 0
void addEvidence(PositionEvidence &posEvidence, const SamRecordView &read,
        int minMappingQuality, int minBaseQuality, BamScanStats &stats)
{
    const int baseQuality = getBaseQuality(posEvidence.position, read);

    bool filtered = false;

    if(read.mappingQuality() < minMappingQuality){
        ++posEvidence.mapQualFiltered;
        filtered = true;
    }

    if(baseQuality < minBaseQuality){
        ++posEvidence.baseQualFiltered;
        filtered = true;
    }

    const int bafPosition = posEvidence.position;
    const int readIndex = read.readPositionAtReferencePosition(bafPosition) - 1; // is 1-based

    // isUltima() 的分支：SEQUENCING_TYPE 非 ULTIMA 時 seqTechFiltered 恆為 0，此處不實作

    ++posEvidence.readDepth;

    if(filtered){
        return;
    }

    if(readIndex >= 0){
        if(!isIndel(bafPosition, readIndex + 1, read)){
            const char baseChar = read.baseAt(readIndex + 1);

            // Java 端在此以 AmberBase.valueOf 解析，遇到非 ACGTN 的 IUPAC 碼會丟例外。
            // C++ 不複製崩潰，改為計數，讓這個已知的不對稱點變成可觀察量。
            if(!isAcgtn(baseChar)){
                ++stats.nonAcgtnBases;
            }

            if(baseChar == posEvidence.ref){
                ++posEvidence.refSupport;
            }else if(baseChar == posEvidence.alt){
                ++posEvidence.altSupport;
            }
        }else{
            ++posEvidence.indelCount;
        }
    }
}

// 對應 RegionTask.processRecord（RegionTask.java:37-61）。
// currentIndex 單調前進，只在 read 起點越過該位點時前進。
void processRecord(RegionTask &task, const SamRecordView &read,
        int minMappingQuality, int minBaseQuality, BamScanStats &stats)
{
    const int alignmentStart = read.alignmentStart();
    const int alignmentEnd = read.alignmentEnd();

    for(std::size_t index = task.currentIndex; index < task.positions.size(); ++index){
        PositionEvidence &posEvidence = *task.positions[index];

        if(alignmentStart > posEvidence.position){
            ++task.currentIndex;
            continue;
        }

        if(alignmentEnd < posEvidence.position){
            break;
        }

        addEvidence(posEvidence, read, minMappingQuality, minBaseQuality, stats);
    }

    // Java 端在此設定 mComplete 並觸發 haltProcessing。behaviour-contract.md §2.4 已論證
    // 該條件在本設定下不可能成立（要觸發它的 read 不會被查詢送進來），故不實作。
}

}

// 共用掃描層（EXP-I02）需要的兩個入口。兩者都只是把既有的匿名 namespace 實作轉出去，
// **一行邏輯都沒有改**——processBam 與 ContigSink 因此呼叫的是同一份程式碼。
bool passesSlicerFilters(const bam1_t *record)
{
    return passesFilters(record);
}

void processRecordForRegion(RegionTask &task, const SamRecordView &read,
        int minMappingQuality, int minBaseQuality, BamScanStats &stats)
{
    processRecord(task, read, minMappingQuality, minBaseQuality, stats);
}

std::vector<RegionTask> populateTaskQueue(
        const std::vector<std::pair<std::string, std::vector<PositionEvidence *>>> &chrPositions,
        int minGap)
{
    std::vector<RegionTask> tasks;

    for(const auto &entry : chrPositions){
        const std::string &chromosome = entry.first;
        const std::vector<PositionEvidence *> &positions = entry.second;

        if(positions.empty()){
            continue;
        }

        tasks.push_back(RegionTask{chromosome, positions[0]->position, positions[0]->position,
                {positions[0]}, 0});

        for(std::size_t i = 1; i < positions.size(); ++i){
            PositionEvidence *posEvidence = positions[i];

            if(tasks.back().end + minGap < posEvidence->position){
                // start a new region
                tasks.push_back(RegionTask{chromosome, posEvidence->position, posEvidence->position,
                        {posEvidence}, 0});
            }else{
                tasks.back().positions.push_back(posEvidence);
                tasks.back().end = std::max(tasks.back().end, posEvidence->position);
            }
        }
    }

    return tasks;
}

BamScanStats processBam(
        const std::string &bamFile, std::vector<RegionTask> &tasks,
        int minMappingQuality, int minBaseQuality, int threads)
{
    const int workerCount = std::max(1, threads);

    // 每個 worker 自帶 samFile / index / bam1_t（htslib 的這些物件非執行緒安全）。
    // 工作分配用一個原子計數器，等同 Java 端的 TaskQueue：
    // 哪個執行緒拿到哪個 region 不影響結果，因為 region 之間的位點互不重疊。
    std::atomic<std::size_t> nextTask{0};
    std::vector<BamScanStats> perWorker(workerCount);
    std::vector<std::string> errors(workerCount);

    const auto work = [&](int workerIndex){
        samFile *in = sam_open(bamFile.c_str(), "r");
        if(in == nullptr){
            errors[workerIndex] = "unable to open bam: " + bamFile;
            return;
        }

        bam_hdr_t *header = sam_hdr_read(in);
        hts_idx_t *index = header == nullptr ? nullptr : sam_index_load(in, bamFile.c_str());
        if(header == nullptr || index == nullptr){
            errors[workerIndex] = "unable to read bam header or index: " + bamFile;
            if(header != nullptr){ sam_hdr_destroy(header); }
            sam_close(in);
            return;
        }

        bam1_t *record = bam_init1();
        BamScanStats &stats = perWorker[workerIndex];

        while(true){
            const std::size_t taskIndex = nextTask.fetch_add(1);
            if(taskIndex >= tasks.size()){
                break;
            }

            RegionTask &task = tasks[taskIndex];

            const int tid = bam_name2id(header, task.chromosome.c_str());
            if(tid < 0){
                continue;
            }

            hts_itr_t *iter = sam_itr_queryi(index, tid, task.start - 1, task.end);
            if(iter == nullptr){
                continue;
            }

            int ret;
            while((ret = sam_itr_next(in, iter, record)) >= 0){
                if(!passesFilters(record)){
                    continue;
                }

                ++stats.recordsConsumed;

                const SamRecordView read(record);
                processRecord(task, read, minMappingQuality, minBaseQuality, stats);
            }

            hts_itr_destroy(iter);

            if(ret < -1){
                errors[workerIndex] = "error reading bam region " + task.chromosome;
                break;
            }
        }

        bam_destroy1(record);
        hts_idx_destroy(index);
        sam_hdr_destroy(header);
        sam_close(in);
    };

    if(workerCount == 1){
        work(0);
    }else{
        std::vector<std::thread> workers;
        workers.reserve(workerCount);
        for(int i = 0; i < workerCount; ++i){
            workers.emplace_back(work, i);
        }
        for(std::thread &worker : workers){
            worker.join();
        }
    }

    for(const std::string &error : errors){
        if(!error.empty()){
            throw std::runtime_error(error);
        }
    }

    // 兩個計數器都是整數加總，與相加順序無關
    BamScanStats stats;
    for(const BamScanStats &s : perWorker){
        stats.recordsConsumed += s.recordsConsumed;
        stats.nonAcgtnBases += s.nonAcgtnBases;
    }

    return stats;
}

}
