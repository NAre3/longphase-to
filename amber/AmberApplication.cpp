// AMBER 4.3 tumor-only 移植：
//   第一段 CP-A1（site 載入）、CP-A2（tumor-only blacklist 過濾）
//   第二段 CP-A2b（region 切分）、CP-A3（per-locus evidence 收集）
//
// 對照對象是 hmftools tag amber-v4.3 的 AmberApplication.loadAmberSites()、
// hetLociTumorOnly()、TumorAnalysis.tumorBAFAndContamination() 與 BamEvidenceReader。
// 本階段刻意不共用 LongPhase-TO 既有的 BAM 掃描，也不接進主程式的建置，理由見研究規格 §0。
//
// 逐一對齊的行為出處：RUN-004/behaviour-contract.md

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "AmberSitesFile.h"
#include "BamEvidenceReader.h"
#include "ChrBaseRegion.h"
#include "CpDump.h"
#include "HumanChromosome.h"
#include "PositionEvidence.h"

namespace {

// 對應 AmberUtils.loadBedFromResource：BED 為 0-based 半開區間，
// 載入時 start+1 轉成 1-based 含端點，end 不變。
std::vector<amber::ChrBaseRegion> loadBed(const std::string &path){
    std::ifstream in(path);
    if(!in){
        throw std::runtime_error("unable to open bed: " + path);
    }

    std::vector<amber::ChrBaseRegion> regions;
    std::string line;

    while(std::getline(in, line)){
        if(line.empty()){
            continue;
        }

        std::istringstream stream(line);
        std::string chromosome;
        int start = 0;
        int end = 0;
        stream >> chromosome >> start >> end;

        regions.push_back({chromosome, start + 1, end});
    }

    return regions;
}

std::string usage(){
    return "usage: amber_port -loci <AmberGermlineSites.tsv.gz> "
           "-tumor_only_excluded_bed <tumorOnlyExcludedSnp.38.bed> -cpdump_dir <dir>\n"
           "                  [-tumor_bam <bam>] [-min_base_quality N] [-min_map_quality N]\n"
           "                  [-debug_only_chr <chr>]\n"
           "\n"
           "-tumor_bam 未給定時只跑到 CP-A2 為止。\n"
           "-debug_only_chr 僅為迭代時縮短週期用的 harness 便利旗標，**不是移植的行為**\n"
           "  （AMBER 的 -specific_chr 並不會限制 loci）。正式記錄的執行必須不帶此旗標。\n";
}

}

int main(int argc, char **argv){
    std::string lociPath;
    std::string bedPath;
    std::string cpDumpDir;
    std::string tumorBam;
    std::string debugOnlyChr;
    int minBaseQuality = amber::DEFAULT_MIN_BASE_QUALITY;
    int minMappingQuality = amber::DEFAULT_MIN_MAPPING_QUALITY;

    for(int i = 1; i < argc - 1; ++i){
        const std::string arg = argv[i];
        if(arg == "-loci"){
            lociPath = argv[++i];
        }else if(arg == "-tumor_only_excluded_bed"){
            bedPath = argv[++i];
        }else if(arg == "-cpdump_dir"){
            cpDumpDir = argv[++i];
        }else if(arg == "-tumor_bam"){
            tumorBam = argv[++i];
        }else if(arg == "-min_base_quality"){
            minBaseQuality = std::stoi(argv[++i]);
        }else if(arg == "-min_map_quality"){
            minMappingQuality = std::stoi(argv[++i]);
        }else if(arg == "-debug_only_chr"){
            debugOnlyChr = argv[++i];
        }
    }

    if(lociPath.empty() || bedPath.empty()){
        std::cerr << usage();
        return 1;
    }

    amber::CpDump::setDir(cpDumpDir);

    // ---- CP-A1：site 載入 ----
    const std::vector<amber::AmberSite> sites = amber::loadAmberSites(lociPath);
    std::fprintf(stderr, "loaded %zu Amber germline sites from %s\n", sites.size(), lociPath.c_str());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(sites.size());
        for(const amber::AmberSite &site : sites){
            rows.push_back({site.chromosome, site.position,
                    site.chromosome + "\t" + std::to_string(site.position) + "\t" + site.ref + "\t" + site.alt});
        }
        amber::CpDump::writeSorted("CP-A1", "chromosome\tposition\tref\talt", rows);
    }

    // ---- CP-A2：tumor-only blacklist 過濾 ----
    // 對應 hetLociTumorOnly()：對每個 site 掃過全部排除區間，命中即計入 numBlackListed。
    // Java 端用的是逐區間線性掃描並在命中時 break，區間只有 32 個，這裡照做，
    // 不改成區間樹——保真度階段以行為一致為先，效能是第二階段的事。
    const std::vector<amber::ChrBaseRegion> excluded = loadBed(bedPath);

    std::vector<amber::AmberSite> retained;
    retained.reserve(sites.size());
    int numBlackListed = 0;

    for(const amber::AmberSite &site : sites){
        bool blacklisted = false;
        for(const amber::ChrBaseRegion &region : excluded){
            if(region.containsPosition(site.chromosome, site.position)){
                blacklisted = true;
                break;
            }
        }

        if(blacklisted){
            ++numBlackListed;
        }else{
            retained.push_back(site);
        }
    }

    std::fprintf(stderr, "removed %d blacklisted loci, %zu remaining\n", numBlackListed, retained.size());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(retained.size());
        for(const amber::AmberSite &site : retained){
            rows.push_back({site.chromosome, site.position,
                    site.chromosome + "\t" + std::to_string(site.position) + "\t" + site.ref + "\t" + site.alt});
        }
        amber::CpDump::writeSorted("CP-A2", "chromosome\tposition\tref\talt", rows);

        std::vector<amber::CpDump::Row> summary;
        summary.push_back({"", 0, "numBlackListed\t" + std::to_string(numBlackListed)});
        summary.push_back({"", 0, "remaining\t" + std::to_string(retained.size())});
        amber::CpDump::write("CP-A2-summary", "field\tvalue", summary);
    }

    if(tumorBam.empty()){
        return 0;
    }

    // ---- CP-A2b / CP-A3：region 切分與 BAM evidence 收集 ----
    //
    // 對應 TumorAnalysis.tumorBAFAndContamination（TumorAnalysis.java:45-98）在 tumor-only 下的
    // 情形：germlineHomLoci 為空（AmberApplication.java:259 傳入空的 multimap），因此
    // chrPositionEvidence 就是 CP-A2 留下的 site 各自建一個歸零的 PositionEvidence
    // （TumorBAF.fromNormal → new PositionEvidence，七個計數器皆 0；TumorBAF.java:34-41）。

    std::vector<amber::PositionEvidence> evidence;
    evidence.reserve(retained.size());

    for(const amber::AmberSite &site : retained){
        if(!debugOnlyChr.empty() && site.chromosome != debugOnlyChr){
            continue;
        }

        amber::PositionEvidence pe;
        pe.chromosome = site.chromosome;
        pe.position = site.position;
        pe.ref = site.ref.empty() ? 'N' : site.ref[0];
        pe.alt = site.alt.empty() ? 'N' : site.alt[0];
        evidence.push_back(std::move(pe));
    }

    // 依染色體分組。Java 端的 key 是 HumanChromosome（去 chr 前綴），region 的染色體字串則是
    // RefGenVersion.versionedChromosome(key)，V38 下為 "chr" + 去前綴名。
    // 分組之間的順序不影響結果（各染色體獨立，且 dump 會再排序）。
    std::vector<std::pair<std::string, std::vector<amber::PositionEvidence *>>> chrPositions;
    std::unordered_map<std::string, std::size_t> chrIndex;

    for(amber::PositionEvidence &pe : evidence){
        const std::string key = amber::stripChrPrefix(pe.chromosome);
        auto found = chrIndex.find(key);

        if(found == chrIndex.end()){
            chrIndex.emplace(key, chrPositions.size());
            chrPositions.emplace_back("chr" + key, std::vector<amber::PositionEvidence *>{&pe});
        }else{
            chrPositions[found->second].second.push_back(&pe);
        }
    }

    // 對應 TumorAnalysis.java:92-95 的 Collections.sort，**在切分之前**執行。
    // 同一組內染色體字串相同，ContigComparator 因此退化為單純比 position；
    // Java 的 List.sort 為 TimSort（穩定），故用 stable_sort。
    for(auto &entry : chrPositions){
        std::stable_sort(entry.second.begin(), entry.second.end(),
                [](const amber::PositionEvidence *a, const amber::PositionEvidence *b){
                    return a->position < b->position;
                });
    }

    std::vector<amber::RegionTask> tasks = amber::populateTaskQueue(chrPositions, amber::BAM_MIN_GAP_START);

    std::fprintf(stderr, "split %zu sites across %zu regions, minGap(%d)\n",
            evidence.size(), tasks.size(), amber::BAM_MIN_GAP_START);

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(tasks.size());
        for(const amber::RegionTask &task : tasks){
            rows.push_back({task.chromosome, task.start,
                    task.chromosome + "\t" + std::to_string(task.start) + "\t"
                            + std::to_string(task.end) + "\t"
                            + std::to_string(task.positions.size())});
        }
        amber::CpDump::writeSorted("CP-A2b", "chromosome\tstart\tend\tpositionCount", rows);
    }

    const amber::BamScanStats stats =
            amber::processBam(tumorBam, tasks, minMappingQuality, minBaseQuality);

    std::fprintf(stderr, "consumed %llu reads, non-ACGTN bases at evaluated positions: %llu\n",
            static_cast<unsigned long long>(stats.recordsConsumed),
            static_cast<unsigned long long>(stats.nonAcgtnBases));

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(evidence.size());
        for(const amber::PositionEvidence &pe : evidence){
            std::string line = pe.chromosome;
            line += "\t" + std::to_string(pe.position);
            line += "\t" + std::string(1, pe.ref);
            line += "\t" + std::string(1, pe.alt);
            line += "\t" + std::to_string(pe.readDepth);
            line += "\t" + std::to_string(pe.indelCount);
            line += "\t" + std::to_string(pe.refSupport);
            line += "\t" + std::to_string(pe.altSupport);
            line += "\t" + std::to_string(pe.baseQualFiltered);
            line += "\t" + std::to_string(pe.mapQualFiltered);
            line += "\t" + std::to_string(pe.seqTechFiltered);
            rows.push_back({pe.chromosome, pe.position, std::move(line)});
        }
        amber::CpDump::writeSorted("CP-A3",
                "chromosome\tposition\tref\talt\treadDepth\tindelCount"
                "\trefSupport\taltSupport\tbaseQualFiltered\tmapQualFiltered\tseqTechFiltered",
                rows);
    }

    return 0;
}
