// COBALT tumor-only whole-genome 的 C++ 移植（purple-port-cobalt-fidelity-v2）。
// 本檔目前涵蓋 EXP-C003（CP-C1..C5）、EXP-C004（CP-C6）與 EXP-C005（CP-C7/C8）。
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

std::string arg(int argc, char **argv, const char *flag, const std::string &fallback = std::string())
{
    for(int i = 1; i + 1 < argc; ++i)
    {
        if(std::strcmp(argv[i], flag) == 0){ return argv[i + 1]; }
    }
    return fallback;
}

CpDump::Row line(const std::string &text){ return CpDump::Row{std::string(), 0, text}; }

}

int main(int argc, char **argv)
{
    const std::string bamPath      = arg(argc, argv, "-tumor_bam");
    const std::string gcProfile    = arg(argc, argv, "-gc_profile");
    const std::string diploidBed   = arg(argc, argv, "-tumor_only_diploid_bed");
    const std::string excludedPath = arg(argc, argv, "-excluded_regions");
    const int threads = std::atoi(arg(argc, argv, "-threads", "1").c_str());
    const int minMappingQuality = std::atoi(arg(argc, argv, "-min_quality", "10").c_str());
    bool includeDuplicates = false;
    for(int i = 1; i < argc; ++i){ if(std::strcmp(argv[i], "-include_duplicates") == 0){ includeDuplicates = true; } }
    if(bamPath.empty() || gcProfile.empty() || excludedPath.empty())
    {
        std::fprintf(stderr,
            "usage: cobalt_port -tumor <id> -tumor_bam <bam> -gc_profile <cnp> "
            "-excluded_regions <tsv> [-tumor_only_diploid_bed <bed.gz>] [-threads N] [-output_dir DIR]\n"
            "  checkpoint dump 由環境變數 COBALT_CPDUMP_DIR 啟用（與 Java 端相同）\n");
        return 2;
    }

    const char *dumpDir = std::getenv("COBALT_CPDUMP_DIR");
    if(dumpDir != nullptr && dumpDir[0] != '\0'){ CpDump::setDir(dumpDir); }

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

    // ---------------- CP-C6 ----------------
    const std::vector<cobalt::DepthReading> depths =
        cobalt::calculateReadDepths(bamPath, chromosomes, partitions, minMappingQuality, includeDuplicates, threads);
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

    std::fprintf(stderr, "cobalt_port: gcReplaced=%ld refLookupFailures=%ld bucketedReadings=%ld\n",
                 buildStats.gcReplacedCount, buildStats.referenceLookupFailures, buildStats.bucketedReadings);
    std::fprintf(stderr, "cobalt_port: depthReadings=%zu threads=%d\n", depths.size(), threads);
    std::fprintf(stderr, "cobalt_port: chromosomes=%zu partitions=%zu gcWindows=%zu diploidEntries=%zu\n",
                 chromosomes.size(), partitions.size(),
                 [&]{ std::size_t n = 0; for(const auto &v : gcData.byChromosome){ n += v.size(); } return n; }(),
                 diploid.totalEntries());
    return 0;
}
