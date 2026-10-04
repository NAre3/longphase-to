#include "ReadDepth.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <thread>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "CobaltConstants.h"

namespace cobalt {

namespace {

// htsjdk AlignmentBlock：只有 M/=/X 產生 block；readBase 由 1 起算、refBase 由 alignmentStart 起算
struct AlignmentBlock
{
    int readStart;        // 1-based，索引到 SEQ（含 soft clip）
    int referenceStart;   // 1-based
    int length;
};

std::vector<AlignmentBlock> alignmentBlocks(const bam1_t *b)
{
    std::vector<AlignmentBlock> blocks;
    const uint32_t *cigar = bam_get_cigar(b);
    int readBase = 1;
    int refBase = static_cast<int>(b->core.pos) + 1;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i)
    {
        const int op = bam_cigar_op(cigar[i]);
        const int len = static_cast<int>(bam_cigar_oplen(cigar[i]));
        switch(op)
        {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                blocks.push_back(AlignmentBlock{readBase, refBase, len});
                readBase += len;
                refBase += len;
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                readBase += len;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                refBase += len;
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
                break;                                  // 兩者皆不前進
            default:
                break;
        }
    }
    return blocks;
}

// BamSlicerFilter.passesFilters（位元組碼逐條，見 behaviour-contract.md §1.2）
inline bool passesFilters(const bam1_t *b, int minMappingQuality)
{
    if(b->core.qual < static_cast<unsigned>(minMappingQuality)){ return false; }   // 嚴格 <
    if((b->core.flag & BAM_FUNMAP) != 0){ return false; }                          // mKeepUnmapped = false
    if((b->core.flag & BAM_FSECONDARY) != 0){ return false; }                      // mKeepSecondaries = false
    if((b->core.flag & BAM_FSUPPLEMENTARY) != 0){ return false; }                  // mKeepSupplementaries = false
    if((b->core.flag & BAM_FDUP) != 0){ return false; }                            // mKeepDuplicates = IncludeDuplicates
    return true;
}

std::string readBasesOf(const bam1_t *b)
{
    const uint8_t *seq = bam_get_seq(b);
    const int n = b->core.l_qseq;
    std::string out(static_cast<std::size_t>(n), 'N');
    for(int i = 0; i < n; ++i){ out[static_cast<std::size_t>(i)] = seq_nt16_str[bam_seqi(seq, i)]; }
    return out;
}

}

void ReadDepthAccumulator::addChromosome(const std::string &chromosome, int chromosomeLength)
{
    if(find(chromosome) != nullptr)
    {
        throw std::runtime_error("duplicate chromosome in accumulator: " + chromosome);   // Validate.isTrue
    }
    const int numWindows = chromosomeLength / mWindowSize;                                 // 整數除法（G10）
    auto c = std::make_unique<ChromosomeWindowCounts>();
    c->numWindows = numWindows;
    c->baseCounts = std::make_unique<std::atomic<int>[]>(static_cast<std::size_t>(numWindows));
    c->gcCounts = std::make_unique<std::atomic<int>[]>(static_cast<std::size_t>(numWindows));
    for(int i = 0; i < numWindows; ++i)
    {
        c->baseCounts[static_cast<std::size_t>(i)].store(0, std::memory_order_relaxed);
        c->gcCounts[static_cast<std::size_t>(i)].store(0, std::memory_order_relaxed);
    }
    mNames.push_back(chromosome);
    mCounts.push_back(std::move(c));
}

const ReadDepthAccumulator::ChromosomeWindowCounts *
ReadDepthAccumulator::find(const std::string &chromosome) const
{
    for(std::size_t i = 0; i < mNames.size(); ++i)
    {
        if(mNames[i] == chromosome){ return mCounts[i].get(); }
    }
    return nullptr;
}

ReadDepthAccumulator::ChromosomeWindowCounts *
ReadDepthAccumulator::find(const std::string &chromosome)
{
    return const_cast<ChromosomeWindowCounts *>(
        static_cast<const ReadDepthAccumulator *>(this)->find(chromosome));
}

void ReadDepthAccumulator::addReadAlignmentToCounts(const std::string &chromosome, int genomeStart,
                                                    int alignmentLength, const std::string &readBases,
                                                    int readStartIndex)
{
    ChromosomeWindowCounts *counts = find(chromosome);
    if(counts == nullptr){ return; }                       // not a chromosome we keep track of

    for(int windowIndex = getWindowIndex(genomeStart); ; ++windowIndex)
    {
        const int windowStart = getGenomePosition(windowIndex);

        if(windowStart >= genomeStart + alignmentLength){ break; }    // 順序：先判這一條
        if(windowIndex >= counts->numWindows){ break; }               // 再判染色體尾端（G10）

        const int startOffset = std::max(0, windowStart - genomeStart);
        const int endOffset = std::max(0, genomeStart + alignmentLength - windowStart - mWindowSize);
        const int numBasesInWindow = alignmentLength - startOffset - endOffset;

        counts->baseCounts[static_cast<std::size_t>(windowIndex)]
            .fetch_add(numBasesInWindow, std::memory_order_relaxed);

        int numGcs = 0;
        for(int i = readStartIndex + startOffset; i < readStartIndex + startOffset + numBasesInWindow; ++i)
        {
            const char base = readBases[static_cast<std::size_t>(i)];
            if(base == 'G' || base == 'C'){ ++numGcs; }     // SequenceUtil.G / .C，只認大寫
        }
        counts->gcCounts[static_cast<std::size_t>(windowIndex)]
            .fetch_add(numGcs, std::memory_order_relaxed);
    }
}

std::vector<DepthReading> ReadDepthAccumulator::getChromosomeReadDepths(const std::string &chromosome) const
{
    const ChromosomeWindowCounts *counts = find(chromosome);
    std::vector<DepthReading> out;
    if(counts == nullptr){ return out; }

    out.reserve(static_cast<std::size_t>(counts->numWindows));
    for(int windowIndex = 0; windowIndex < counts->numWindows; ++windowIndex)
    {
        const double basesCount = counts->baseCounts[static_cast<std::size_t>(windowIndex)]
                                      .load(std::memory_order_relaxed);
        const double depth = basesCount / mWindowSize;
        // basesCount == 0 時為 NaN；保留「先算出 NaN 再被覆蓋」的形式（見對照表 §5）
        const double gcContent = counts->gcCounts[static_cast<std::size_t>(windowIndex)]
                                     .load(std::memory_order_relaxed) / basesCount;
        DepthReading r;
        r.chromosome = chromosome;
        r.startPosition = getGenomePosition(windowIndex);
        r.readDepth = depth;
        r.readGcContent = (depth == 0) ? 0 : gcContent;    // DepthReading 建構子（G11）
        out.push_back(std::move(r));
    }
    return out;
}

// 共用掃描層（EXP-I03）。函式體就是 calculateReadDepths 內層 while 迴圈的那一段，
// **逐字相同**，只是把 region 固定成 {chromosome, 1, contigLength}（design.md E1）。
// 兩條路徑因此走同一份過濾與裁切邏輯。
DepthSink::DepthSink(ReadDepthAccumulator &accumulator, std::string chromosome, int contigLength,
                     int minMappingQuality, bool includeDuplicates)
    : mAccumulator(accumulator), mChromosome(std::move(chromosome)), mContigLength(contigLength),
      mMinMappingQuality(minMappingQuality), mIncludeDuplicates(includeDuplicates)
{
}

void DepthSink::consume(const bam1_t *rec)
{
    // COBALT 自己的 filter。**不可與 AMBER 共用**：COBALT 排除 qual < 10（嚴格），
    // AMBER 在 read 層完全不過濾（其 MAPQ 是在 addEvidence 內計數而非丟棄）。
    // 把 MAPQ 上提到共用層會讓 AMBER 的 ReadDepth 少算而其餘欄位不變——ALT-1 的表徵。
    if(!passesFilters(rec, mMinMappingQuality)){ return; }

    // processRead 的第二道檢查（與 calculateReadDepths 同）
    if(mIncludeDuplicates){ if(bam_aux_get(rec, "CR") != nullptr){ return; } }
    else                  { if((rec->core.flag & BAM_FDUP) != 0){ return; } }

    const std::string bases = readBasesOf(rec);
    for(const AlignmentBlock &block : alignmentBlocks(rec))
    {
        const int genomeStart = std::max(block.referenceStart, 1);
        const int length = std::min(block.referenceStart + block.length, mContigLength + 1) - genomeStart;
        if(length <= 0){ continue; }
        const int readStartIndex = (block.readStart - 1) + (genomeStart - block.referenceStart);
        mAccumulator.addReadAlignmentToCounts(mChromosome, genomeStart, length, bases, readStartIndex);
    }
}

std::vector<DepthReading> calculateReadDepths(const std::string &bamPath,
                                              const std::vector<ChromosomeSpec> &chromosomes,
                                              const std::vector<lp::ChrBaseRegion> &partitions,
                                              int minMappingQuality,
                                              bool includeDuplicates,
                                              int threads)
{
    ReadDepthAccumulator accumulator(WINDOW_SIZE);
    for(const ChromosomeSpec &c : chromosomes){ accumulator.addChromosome(c.name, c.length); }

    std::atomic<std::size_t> next(0);
    std::atomic<bool> failed(false);
    const int workers = std::max(1, threads);

    auto worker = [&]()
    {
        samFile *fp = sam_open(bamPath.c_str(), "r");
        if(fp == nullptr){ failed.store(true); return; }
        hts_idx_t *idx = sam_index_load(fp, bamPath.c_str());
        sam_hdr_t *hdr = sam_hdr_read(fp);
        if(idx == nullptr || hdr == nullptr){ failed.store(true); sam_close(fp); return; }
        bam1_t *rec = bam_init1();

        while(true)
        {
            const std::size_t i = next.fetch_add(1);
            if(i >= partitions.size()){ break; }
            const lp::ChrBaseRegion &region = partitions[i];
            const int tid = sam_hdr_name2tid(hdr, region.chromosome.c_str());
            if(tid < 0){ continue; }

            // ChrBaseRegion 為 1-based 含端點；htslib 為 0-based 半開
            hts_itr_t *it = sam_itr_queryi(idx, tid, region.start - 1, region.end);
            if(it == nullptr){ continue; }

            while(sam_itr_next(fp, it, rec) >= 0)
            {
                if(!passesFilters(rec, minMappingQuality)){ continue; }
                // processRead 的第二道檢查（本 run includeDuplicates = false，與上一行重複）
                if(includeDuplicates){ if(bam_aux_get(rec, "CR") != nullptr){ continue; } }
                else                 { if((rec->core.flag & BAM_FDUP) != 0){ continue; } }

                const std::string bases = readBasesOf(rec);
                for(const AlignmentBlock &block : alignmentBlocks(rec))
                {
                    const int genomeStart = std::max(block.referenceStart, region.start);
                    const int length = std::min(block.referenceStart + block.length, region.end + 1) - genomeStart;
                    if(length <= 0){ continue; }
                    const int readStartIndex = (block.readStart - 1) + (genomeStart - block.referenceStart);
                    accumulator.addReadAlignmentToCounts(region.chromosome, genomeStart, length,
                                                         bases, readStartIndex);
                }
            }
            hts_itr_destroy(it);
        }

        bam_destroy1(rec);
        sam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp);
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(workers));
    for(int i = 0; i < workers; ++i){ pool.emplace_back(worker); }
    for(std::thread &t : pool){ t.join(); }
    if(failed.load()){ throw std::runtime_error("bam reading failed: " + bamPath); }

    // 輸出順序 = CP-C1 的染色體順序，各染色體內 window 索引昇冪
    std::vector<DepthReading> out;
    for(const ChromosomeSpec &c : chromosomes)
    {
        std::vector<DepthReading> v = accumulator.getChromosomeReadDepths(c.name);
        out.insert(out.end(), std::make_move_iterator(v.begin()), std::make_move_iterator(v.end()));
    }
    return out;
}

}
