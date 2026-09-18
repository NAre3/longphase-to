// COBALT tumor-only whole-genome 的 C++ 移植（purple-port-cobalt-fidelity-v2）。
// 本檔目前涵蓋 EXP-C003（CP-C1..C5）、EXP-C004（CP-C6）、EXP-C005（CP-C7/C8）
// EXP-C006（CP-C9/C10）、EXP-C007（CP-C11/C12 與兩個 stage 輸出）
// 與 EXP-C008（CP-C13/C13b/C13c/C14 與 cobalt.ratio.pcf）。
//
// 行為出處逐條見 runs/RUN-C003/behaviour-contract.md。

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "CobaltConstants.h"
#include "GcProfile.h"
#include "BamRatio.h"
#include "CobaltOutput.h"
#include "Consolidation.h"
#include "Segmentation.h"
#include "CobaltWindow.h"
#include "GcBuckets.h"
#include "ReadDepth.h"
#include "Regions.h"
#include "WindowStatuses.h"
#include "../common/ChrBaseRegion.h"
#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"
#include "../common/Segmentation.h"   // lp::chromosomeOrdinal

using lp::CpDump;

namespace {

// BamReadCounter.loadChromosomes：依 header @SQ 順序，只留 HumanChromosome.contains 為真者
std::vector<cobalt::ChromosomeSpec> loadChromosomes(const std::string &bamPath)
{
    samFile *fp = sam_open(bamPath.c_str(), "r");
    if(fp == nullptr){ std::fprintf(stderr, "cannot open bam: %s\n", bamPath.c_str()); std::exit(2); }
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if(hdr == nullptr){ std::fprintf(stderr, "cannot read bam header\n"); std::exit(2); }

    std::vector<cobalt::ChromosomeSpec> out;
    const int n = sam_hdr_nref(hdr);
    for(int i = 0; i < n; ++i)
    {
        const char *name = sam_hdr_tid2name(hdr, i);
        if(name == nullptr){ continue; }
        const std::string sequenceName(name);
        if(!lp::isHumanChromosome(sequenceName)){ continue; }
        // SpecificChrRegions 未設定，該分支不執行
        out.push_back(cobalt::ChromosomeSpec{sequenceName, static_cast<int>(sam_hdr_tid2len(hdr, i))});
    }
    sam_hdr_destroy(hdr);
    sam_close(fp);
    return out;
}

// BamReadCounter.partitionGenome
std::vector<lp::ChrBaseRegion> partitionGenome(const std::vector<cobalt::ChromosomeSpec> &chromosomes)
{
    std::vector<lp::ChrBaseRegion> partitions;
    for(const cobalt::ChromosomeSpec &c : chromosomes)
    {
        for(int startPos = 1; startPos < c.length; startPos += cobalt::PARTITION_SIZE)
        {
            const int endPos = std::min(startPos + cobalt::PARTITION_SIZE - 1, c.length);
            partitions.push_back(lp::ChrBaseRegion{c.name, startPos, endPos});
        }
    }
    return partitions;
}

CpDump::Row line(const std::string &text){ return CpDump::Row{std::string(), 0, text}; }

}

#include "CobaltPipeline.h"

namespace cobalt {

// 以下兩個函式的內容是從 CobaltApplication.cpp 的 main() **逐字搬過來的**，
// 只做了兩件事：
//   1. 開頭加一段 const 參考別名，讓搬過來的程式碼不必改寫變數名；
//   2. 原本以 arg(argc, argv, ...) 取得的 -output_dir / -tumor / -pcf_gamma
//      改為讀 cfg 的同名欄位（預設值不變）。
// 除此之外沒有任何邏輯變更。行為出處的註解全部原樣保留。

PrescanResult prescan(const PipelineConfig &cfg)
{
    const std::string &bamPath = cfg.bamPath;
    const std::string &gcProfile = cfg.gcProfile;
    const std::string &diploidBed = cfg.diploidBed;
    const std::string &excludedPath = cfg.excludedPath;

    PrescanResult result;

    // ---------------- CP-C1 / CP-C2 ----------------
    const std::vector<cobalt::ChromosomeSpec> chromosomes = loadChromosomes(bamPath);
    const std::vector<lp::ChrBaseRegion> partitions = partitionGenome(chromosomes);

    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        for(std::size_t i = 0; i < chromosomes.size(); ++i)
        {
            rows.push_back(line(std::to_string(i) + "\t" + chromosomes[i].name + "\t"
                                + std::to_string(chromosomes[i].length)));
        }
        CpDump::write("CP-C1", "idx\tchromosome\tlength", rows);

        rows.clear();
        for(std::size_t i = 0; i < partitions.size(); ++i)
        {
            const lp::ChrBaseRegion &r = partitions[i];
            rows.push_back(line(std::to_string(i) + "\t" + r.chromosome + "\t"
                                + std::to_string(r.start) + "\t" + std::to_string(r.end)));
        }
        CpDump::write("CP-C2", "idx\tchromosome\tstart\tend", rows);
    }

    // ---------------- CP-C3 ----------------
    const cobalt::GcProfileData gcData = cobalt::loadGcProfile(gcProfile);
    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        for(std::size_t ci = 0; ci < gcData.chromosomes.size(); ++ci)
        {
            for(const cobalt::GcProfile &g : gcData.byChromosome[ci])
            {
                rows.push_back(line(g.chromosome + "\t" + std::to_string(g.start) + "\t"
                                    + std::to_string(g.end) + "\t" + CpDump::num(g.gcContent) + "\t"
                                    + CpDump::num(g.nonNPercentage) + "\t"
                                    + CpDump::num(g.mappablePercentage) + "\t"
                                    + (g.isMappable() ? "true" : "false")));
            }
        }
        CpDump::write("CP-C3", "chromosome\tstart\tend\tgcContent\tnonNPercentage\tmappablePercentage\tisMappable", rows);
    }

    // ---------------- CP-C4 ----------------
    const std::vector<cobalt::ExcludedRegion> excluded = cobalt::loadExcludedRegions(excludedPath);
    cobalt::DiploidRegions diploid;
    if(!diploidBed.empty()){ diploid = cobalt::loadDiploidRegions(diploidBed); }

    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        for(const cobalt::ExcludedRegion &r : excluded)
        {
            // Java 端行尾多一個 \t（空的第五欄）
            rows.push_back(line("excluded\t" + r.chromosome + "\t" + std::to_string(r.start)
                                + "\t" + std::to_string(r.end) + "\t"));
        }
        for(std::size_t ci = 0; ci < diploid.chromosomes.size(); ++ci)
        {
            const std::vector<cobalt::DiploidStatus> &ds = diploid.byChromosome[ci];
            for(std::size_t i = 0; i < ds.size(); ++i)
            {
                rows.push_back(line("diploid\t" + diploid.chromosomes[ci] + "\t"
                                    + std::to_string(i) + "\t\t" + (ds[i].isDiploid ? "true" : "false")));
            }
        }
        CpDump::write("CP-C4", "kind\tchromosome\ta\tb\tisDiploid", rows);

        std::vector<CpDump::Row> scalars;
        scalars.push_back(line("excludedRegionCount\t" + std::to_string(excluded.size())));
        scalars.push_back(line("diploidChromosomeCount\t" + std::to_string(diploid.chromosomes.size())));
        scalars.push_back(line("diploidEntryCount\t" + std::to_string(diploid.totalEntries())));
        CpDump::write("CP-C4-summary", "field\tvalue", scalars);
    }

    // ---------------- CP-C5 ----------------
    const cobalt::WindowStatuses statuses = cobalt::buildWindowStatuses(gcData, excluded, diploid);
    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        for(std::size_t ci = 0; ci < statuses.chromosomes.size(); ++ci)
        {
            for(const cobalt::WindowStatus &w : statuses.byChromosome[ci])
            {
                rows.push_back(line(w.chromosome + "\t" + std::to_string(w.start) + "\t"
                                    + std::to_string(w.end) + "\t"
                                    + (w.excluded ? "true" : "false") + "\t"
                                    + (w.unmappable ? "true" : "false") + "\t"
                                    + (w.nonDiploid ? "true" : "false") + "\t"
                                    + (w.maskedOut() ? "true" : "false")));
            }
        }
        CpDump::write("CP-C5", "chromosome\tstart\tend\texcluded\tunmappable\tnonDiploid\tmaskedOut", rows);
    }

    result.chromosomes = std::move(const_cast<std::vector<cobalt::ChromosomeSpec>&>(chromosomes));
    result.partitions = std::move(const_cast<std::vector<lp::ChrBaseRegion>&>(partitions));
    result.gcData = std::move(const_cast<cobalt::GcProfileData&>(gcData));
    result.excluded = std::move(const_cast<std::vector<cobalt::ExcludedRegion>&>(excluded));
    result.diploid = std::move(diploid);
    result.statuses = std::move(const_cast<cobalt::WindowStatuses&>(statuses));
    return result;
}

PostscanResult postscan(const PipelineConfig &cfg, const PrescanResult &pre,
        const std::vector<DepthReading> &depths)
{
    const std::vector<cobalt::ChromosomeSpec> &chromosomes = pre.chromosomes;
    const cobalt::GcProfileData &gcData = pre.gcData;
    const cobalt::WindowStatuses &statuses = pre.statuses;
    // 只被最後那幾行摘要 fprintf 用到，但仍必須存在，否則摘要的數字會變
    const std::vector<lp::ChrBaseRegion> &partitions = pre.partitions;
    const cobalt::DiploidRegions &diploid = pre.diploid;
    const int threads = cfg.threads;
    (void)chromosomes; (void)gcData; (void)statuses; (void)threads;
    (void)partitions; (void)diploid;

    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        rows.reserve(depths.size());
        for(const cobalt::DepthReading &d : depths)
        {
            rows.push_back(line(d.chromosome + "\t" + std::to_string(d.startPosition) + "\t"
                                + CpDump::num(d.readDepth) + "\t" + CpDump::num(d.readGcContent)));
        }
        CpDump::write("CP-C6", "chromosome\tposition\treadDepth\treadGcContent", rows);
    }

    // ---------------- CP-C7 / CP-C8 ----------------
    cobalt::GcPailsList pails;
    cobalt::WindowBuildStats buildStats;
    std::vector<cobalt::CobaltWindow> windows =
        cobalt::buildWindows(depths, statuses, gcData, pails, buildStats);
    const cobalt::GcBucketStatistics bucketStats(pails, cobalt::GC_BUCKET_MIN, cobalt::GC_BUCKET_MAX);

    if(CpDump::enabled())
    {
        // dump 端依 HumanChromosome ordinal 排序，染色體內依 position 昇冪
        // （Java 端 BamCalculation.java:63 的 sortedKeys + 插入順序）
        std::vector<const cobalt::CobaltWindow *> ordered;
        ordered.reserve(windows.size());
        for(const cobalt::CobaltWindow &w : windows){ ordered.push_back(&w); }
        std::stable_sort(ordered.begin(), ordered.end(),
                         [](const cobalt::CobaltWindow *a, const cobalt::CobaltWindow *b){
            const int oa = lp::chromosomeOrdinal(a->chromosomeShort);
            const int ob = lp::chromosomeOrdinal(b->chromosomeShort);
            if(oa != ob){ return oa < ob; }
            return a->position < b->position;
        });

        std::vector<CpDump::Row> rows;
        rows.reserve(ordered.size());
        for(const cobalt::CobaltWindow *w : ordered)
        {
            rows.push_back(line(w->chromosomeShort + "\t" + std::to_string(w->position) + "\t"
                                + CpDump::num(w->readDepth) + "\t" + CpDump::num(w->gcContent) + "\t"
                                + std::to_string(w->gcBucket) + "\t"
                                + (w->isInExcludedRegion ? "true" : "false") + "\t"
                                + (w->isInTargetRegion ? "true" : "false") + "\t"
                                + (w->gcReplaced ? "true" : "false")));
        }
        CpDump::write("CP-C7",
            "chromosome\tposition\treadDepth\tgcContent\tgcBucket\tisInExcludedRegion\tisInTargetRegion\tgcReplaced", rows);

        rows.clear();
        for(int i = 0; i < 101; ++i)
        {
            const cobalt::GcPail &pail = pails.buckets()[static_cast<std::size_t>(i)];
            rows.push_back(line(std::to_string(i) + "\t" + std::to_string(pail.readingCount()) + "\t"
                                + CpDump::num(pail.median()) + "\t"
                                + CpDump::num(bucketStats.medianReadDepthForBucket(i))));
        }
        CpDump::write("CP-C8", "bucket\treadingCount\tpailMedian\tmeanDepth", rows);

        std::vector<CpDump::Row> scalars;
        scalars.push_back(line("gcBucketMin\t" + std::to_string(cobalt::GC_BUCKET_MIN)));
        scalars.push_back(line("gcBucketMax\t" + std::to_string(cobalt::GC_BUCKET_MAX)));
        CpDump::write("CP-C8-summary", "field\tvalue", scalars);
    }

    // ---------------- CP-C9 / CP-C10 ----------------
    // 走訪順序：依 CP-C1 的染色體順序、window 昇冪。
    // **這與 Java 的 ListMultimap 鍵走訪順序（identity hash，G28）不同，且刻意不重現**
    // ——見 spec v2 §1.2 DECISION D2 與 behaviour-contract-mean.md §3.4。
    // 後果：readDepthMean 與 Java 相差約 1e-15（相對），落在 §5.1 第二層容差內。
    std::vector<cobalt::BamRatio> ratios;
    ratios.reserve(windows.size());
    {
        std::vector<const cobalt::CobaltWindow *> ordered;
        ordered.reserve(windows.size());
        for(const cobalt::CobaltWindow &w : windows){ ordered.push_back(&w); }
        std::stable_sort(ordered.begin(), ordered.end(),
                         [](const cobalt::CobaltWindow *a, const cobalt::CobaltWindow *b){
            const int oa = lp::chromosomeOrdinal(a->chromosomeShort);
            const int ob = lp::chromosomeOrdinal(b->chromosomeShort);
            if(oa != ob){ return oa < ob; }
            return a->position < b->position;
        });
        for(const cobalt::CobaltWindow *w : ordered)
        {
            cobalt::BamRatio r = cobalt::BamRatio::fromWindow(*w);
            // normaliseForGc：excluded window 的 bucket 為 null -> isAllowed 為假 -> -1（G14）
            const cobalt::GcPail *pail = (w->gcBucket < 0)
                ? nullptr
                : &pails.buckets()[static_cast<std::size_t>(w->gcBucket)];
            r.normaliseForGc(bucketStats.medianReadDepth(pail));
            r.applyEnrichment(1.0);                 // WholeGenome：enrichmentQuotient 恆 1.0（G15）
            ratios.push_back(std::move(r));
        }
    }

    auto dumpRatios = [&](const char *id, bool withStatsFlags)
    {
        std::vector<CpDump::Row> rows;
        rows.reserve(ratios.size());
        for(const cobalt::BamRatio &r : ratios)
        {
            std::string s = r.chromosomeShort + "\t" + std::to_string(r.position) + "\t"
                          + CpDump::num(r.readDepth()) + "\t" + CpDump::num(r.ratio()) + "\t"
                          + CpDump::num(r.gcContent());
            if(withStatsFlags)
            {
                s += std::string("\t") + (r.included() ? "true" : "false") + "\t"
                   + (cobalt::ReadDepthStatisticsNormaliser::includedInStats(r) ? "true" : "false");
            }
            rows.push_back(line(s));
        }
        std::string header = "chromosome\tposition\treadDepth\tratio\tgcContent";
        if(withStatsFlags){ header += "\tincluded\tincludedInStats"; }
        CpDump::write(id, header, rows);
    };

    if(CpDump::enabled()){ dumpRatios("CP-C9", true); }

    cobalt::ReadDepthStatisticsNormaliser meanNormaliser;
    for(const cobalt::BamRatio &r : ratios){ meanNormaliser.recordValue(r); }
    meanNormaliser.dataCollectionFinished();
    for(cobalt::BamRatio &r : ratios){ meanNormaliser.normalise(r); }

    if(CpDump::enabled())
    {
        dumpRatios("CP-C10", false);
        std::vector<CpDump::Row> scalars;
        scalars.push_back(line("readDepthMean\t" + CpDump::num(meanNormaliser.readDepthMean())));
        scalars.push_back(line("readDepthMedian\t" + CpDump::num(meanNormaliser.readDepthMedian())));
        scalars.push_back(line("sampleCount\t" + std::to_string(meanNormaliser.sampleCount())));
        CpDump::write("CP-C10-summary", "field\tvalue", scalars);
    }

    // ---------------- CP-C11 / CP-C12 與 stage 輸出 ----------------
    const cobalt::ConsolidatorChoice consolidator =
        cobalt::chooseConsolidator(meanNormaliser.readDepthMedian());
    if(consolidator.className == "LowCoverageConsolidator")
    {
        // **本分支在 dev 樣本上未被執行**（NoOp）。實作存在是為了讓 EXP-C009 能在
        // Java 對照下驗證它，而非默默算下去——CP-C11-summary 的類別名是第一層零容差欄位。
        cobalt::applyLowCoverageConsolidation(ratios, consolidator.consolidationCount);
    }
    // MegaBaseScaleNormaliser 與 FinalNormaliser 在 tumor-only whole-genome 下皆為
    // DoNothingNormaliser（G15）：兩趟 forEach 皆為空實作，不改變任何值。

    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        rows.reserve(ratios.size());
        for(const cobalt::BamRatio &r : ratios)
        {
            rows.push_back(line(r.chromosomeShort + "\t" + std::to_string(r.position) + "\t"
                                + CpDump::num(r.readDepth()) + "\t" + CpDump::num(r.ratio()) + "\t"
                                + CpDump::num(r.gcContent()) + "\t"
                                + CpDump::num(r.diploidAdjustedRatio())));
        }
        CpDump::write("CP-C11",
            "chromosome\tposition\treadDepth\tratio\tgcContent\tdiploidAdjustedRatio", rows);

        std::vector<CpDump::Row> scalars;
        scalars.push_back(line("medianReadDepth\t" + CpDump::num(meanNormaliser.readDepthMedian())));
        scalars.push_back(line("consolidator\t" + consolidator.className));
        scalars.push_back(line("megaBaseScaleNormaliser\tDoNothingNormaliser"));
        scalars.push_back(line("finalNormaliser\tDoNothingNormaliser"));
        CpDump::write("CP-C11-summary", "field\tvalue", scalars);
    }

    std::vector<cobalt::CobaltRatio> collated = cobalt::collateResults(ratios);

    if(CpDump::enabled())
    {
        std::vector<CpDump::Row> rows;
        rows.reserve(collated.size());
        for(const cobalt::CobaltRatio &c : collated)
        {
            // CP-C12 用**物件欄序**（與檔案欄序不同，見 behaviour-contract.md §5.1）
            rows.push_back(line(c.chromosome + "\t" + std::to_string(c.position) + "\t"
                                + CpDump::num(c.referenceReadDepth) + "\t" + CpDump::num(c.referenceGCRatio) + "\t"
                                + CpDump::num(c.referenceGCContent) + "\t" + CpDump::num(c.referenceGCDiploidRatio) + "\t"
                                + CpDump::num(c.tumorReadDepth) + "\t" + CpDump::num(c.tumorGCRatio) + "\t"
                                + CpDump::num(c.tumorGcContent)));
        }
        CpDump::write("CP-C12", "chromosome\tposition\treferenceReadDepth\treferenceGCRatio\t"
            "referenceGCContent\treferenceGCDiploidRatio\ttumorReadDepth\ttumorGCRatio\ttumorGCContent", rows);
    }

    // ---- stage 輸出 ----
    // 原本無條件寫出；加上 outputDir 與 writeStageOutputs 兩個條件後，
    // 既有呼叫端（outputDir 預設 "."、writeStageOutputs 預設 true）行為不變。
    if(!cfg.outputDir.empty() && cfg.writeStageOutputs){
        const std::string &outDir = cfg.outputDir;
        const std::string &tumorId = cfg.sampleId;
        cobalt::writeCobaltRatioFile(outDir + "/" + tumorId + ".cobalt.ratio.tsv.gz", collated);
        cobalt::writeGcMedianFile(outDir + "/" + tumorId + ".cobalt.gc.median.tsv",
                                  meanNormaliser.readDepthMean(), meanNormaliser.readDepthMedian(), bucketStats);
    }

    // ---------------- CP-C13 / C13b / C13c / C14 與 cobalt.ratio.pcf ----------------
    // C++ 端的分段是**單執行緒**（Java 的 PerArmSegmenter.getSegmentation 吃 executor，
    // C++ 未平行化）。依票面，本票的執行緒掃描結果僅涵蓋單執行緒分段，須在 manifest 明寫。
    const double pcfGamma = cfg.pcfGamma;
    cobalt::SegmentationResult seg = cobalt::segmentRatios(collated, pcfGamma);

    auto armLabel = [](const cobalt::ArmData &a){
        return "ChrArm[chromosome=" + a.chromosomeShort + ", arm=" + std::string(1, a.arm) + "]";
    };

    // **dump 的 arm 順序與 seg.arms 的順序不同。**
    // seg.arms 依 ChrArm.compareTo（染色體 enum 序）排序，這是 SegmentsFile.write 走的順序，
    // 也是 cobalt.ratio.pcf 的列順序。
    // 但 Java 的 dump 另外以 `arms.sort(Comparator.comparing(ChrArm::toString))` 排——
    // 那是**對整個標籤字串做字典序**，故為 1P 1Q 10P 10Q 11P … 2P 2Q …
    // 兩者不同；evaluator 按記錄鍵 join 故判定不受影響，但對齊之後 dump 才能逐位元組比對。
    std::vector<const cobalt::ArmData *> dumpOrder;
    dumpOrder.reserve(seg.arms.size());
    for(const cobalt::ArmData &a : seg.arms){ dumpOrder.push_back(&a); }
    std::stable_sort(dumpOrder.begin(), dumpOrder.end(),
                     [&](const cobalt::ArmData *a, const cobalt::ArmData *b){
        return armLabel(*a) < armLabel(*b);
    });

    if(CpDump::enabled())
    {
        // ---- CP-C13：valueForSegmentation、rawValue、floored ----
        std::vector<CpDump::Row> rows;
        for(const cobalt::ArmData *ap : dumpOrder)
        {
            const cobalt::ArmData &a = *ap;
            const std::string label = armLabel(a);
            for(std::size_t i = 0; i < a.valuesForSegmentation.size(); ++i)
            {
                rows.push_back(line(label + "\t" + std::to_string(i) + "\t"
                    + CpDump::num(a.valuesForSegmentation[i]) + "\t"
                    + CpDump::num(a.rawValues[i]) + "\t"
                    + (a.rawValues[i] < 0.001 ? "true" : "false")));
            }
        }
        CpDump::write("CP-C13", "chrArm\tidx\tvalueForSegmentation\trawValue\tfloored", rows);

        std::vector<CpDump::Row> scalars;
        scalars.push_back(line("totalCount\t" + std::to_string(seg.totalCount)));
        scalars.push_back(line("uniformPenaltyThreshold\t" + std::to_string(seg.uniformPenaltyThreshold)));
        scalars.push_back(line("penaltyMode\t" + seg.penaltyMode));
        scalars.push_back(line("gamma\t" + CpDump::num(seg.gamma)));
        scalars.push_back(line(std::string("isWindowed\t") + (seg.isWindowed ? "true" : "false")));
        scalars.push_back(line("armCount\t" + std::to_string(seg.arms.size())));
        CpDump::write("CP-C13-summary", "field\tvalue", scalars);

        // ---- CP-C13b：Gamma 的中間量（長格式）----
        rows.clear();
        for(const cobalt::ArmData *ap : dumpOrder)
        {
            const cobalt::ArmData &a = *ap;
            const std::string label = armLabel(a);
            rows.push_back(line(label + "\tn\t" + CpDump::num(static_cast<double>(a.trace.n))));
            rows.push_back(line(label + "\tfilterWidth\t" + CpDump::num(static_cast<double>(a.trace.filterWidth))));
            rows.push_back(line(label + "\tmad\t" + CpDump::num(a.trace.mad)));
            rows.push_back(line(label + "\tpenalty\t" + CpDump::num(a.trace.penalty)));
            rows.push_back(line(label + "\tnormalised\t" + CpDump::num(a.trace.normalised ? 1.0 : 0.0)));
            rows.push_back(line(label + "\tsegmentPenaltyUsed\t" + CpDump::num(a.trace.penalty)));
            for(std::size_t i = 0; i < a.trace.runningMedians.size(); ++i)
            {
                rows.push_back(line(label + "\trunmed[" + std::to_string(i) + "]\t"
                    + CpDump::num(a.trace.runningMedians[i])));
            }
        }
        CpDump::write("CP-C13b", "chrArm\tfield\tvalue", rows);

        // ---- CP-C13c：PiecewiseConstantFit 自己的 lengths/startPositions/means ----
        rows.clear();
        for(const cobalt::ArmData *ap : dumpOrder)
        {
            const cobalt::ArmData &a = *ap;
            const std::string label = armLabel(a);
            for(std::size_t i = 0; i < a.fit.lengths.size(); ++i)
            {
                rows.push_back(line(label + "\t" + std::to_string(i) + "\t"
                    + std::to_string(a.fit.lengths[i]) + "\t"
                    + std::to_string(a.fit.startPositions[i]) + "\t"
                    + CpDump::num(a.pcfMeans[i])));
            }
        }
        CpDump::write("CP-C13c", "chrArm\tidx\tlength\tstartPosition\tpcfMean", rows);

        // ---- CP-C14：最終 segment ----
        rows.clear();
        for(const cobalt::ArmData *ap : dumpOrder)
        {
            const cobalt::ArmData &a = *ap;
            const std::string label = armLabel(a);
            for(std::size_t i = 0; i < a.segments.size(); ++i)
            {
                const cobalt::PcfSegmentOut &sg = a.segments[i];
                rows.push_back(line(label + "\t" + std::to_string(i) + "\t" + sg.chromosome + "\t"
                    + std::to_string(sg.start) + "\t" + std::to_string(sg.end) + "\t"
                    + CpDump::num(sg.meanRatio)));
            }
        }
        CpDump::write("CP-C14", "chrArm\tidx\tchromosome\tstart\tend\tmeanRatio", rows);
    }

    if(!cfg.outputDir.empty() && cfg.writeStageOutputs){
        const std::string &outDir = cfg.outputDir;
        const std::string &tumorId = cfg.sampleId;
        cobalt::writeSegmentsFile(outDir + "/" + tumorId + ".cobalt.ratio.pcf", seg);
    }

    std::fprintf(stderr, "cobalt_port: segmentation arms=%zu totalCount=%d penaltyMode=%s gamma=%g\n",
                 seg.arms.size(), seg.totalCount, seg.penaltyMode.c_str(), seg.gamma);
    std::fprintf(stderr, "cobalt_port: consolidator=%s count=%d\n",
                 consolidator.className.c_str(), consolidator.consolidationCount);
    std::fprintf(stderr, "cobalt_port: readDepthMean=%.17g readDepthMedian=%.17g sampleCount=%ld\n",
                 meanNormaliser.readDepthMean(), meanNormaliser.readDepthMedian(), meanNormaliser.sampleCount());
    std::fprintf(stderr, "cobalt_port: gcReplaced=%ld refLookupFailures=%ld bucketedReadings=%ld\n",
                 buildStats.gcReplacedCount, buildStats.referenceLookupFailures, buildStats.bucketedReadings);
    std::fprintf(stderr, "cobalt_port: depthReadings=%zu threads=%d\n", depths.size(), threads);
    std::fprintf(stderr, "cobalt_port: chromosomes=%zu partitions=%zu gcWindows=%zu diploidEntries=%zu\n",
                 chromosomes.size(), partitions.size(),
                 [&]{ std::size_t n = 0; for(const auto &v : gcData.byChromosome){ n += v.size(); } return n; }(),
                 diploid.totalEntries());

    PostscanResult result;
    result.ratios = std::move(collated);
    result.segmentation = std::move(seg);
    return result;
}

}
