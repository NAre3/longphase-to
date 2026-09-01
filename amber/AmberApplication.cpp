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
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "AmberOutput.h"
#include "AmberSitesFile.h"
#include "BamEvidenceReader.h"
#include "ChrBaseRegion.h"
#include "CpDump.h"
#include "HumanChromosome.h"
#include "NoiseFloor.h"
#include "PositionEvidence.h"
#include "Segmentation.h"
#include "TumorFilters.h"

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
           "                  [-output_dir <dir>] [-tumor <sampleId>]\n"
           "                  [-debug_only_chr <chr>]\n"
           "\n"
           "-tumor_bam 未給定時只跑到 CP-A2 為止。\n"
           "-output_dir 與 -tumor 同時給定時寫出 amber.baf.tsv.gz 與 amber.qc。\n"
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
    std::string outputDir;
    std::string sampleId;
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
        }else if(arg == "-output_dir"){
            outputDir = argv[++i];
        }else if(arg == "-tumor"){
            sampleId = argv[++i];
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

    // ---- CP-A4：IndelCount == 0 的保留集合 ----
    //
    // 對應 TumorAnalysis.java:100-102：tumorBAFs 中 IndelCount 為 0 者放進 mBafs。
    // 注意 CP-A3 是**全部**位點，CP-A4 才是第一次縮減。

    std::vector<const amber::PositionEvidence *> retainedEvidence;
    retainedEvidence.reserve(evidence.size());

    for(const amber::PositionEvidence &pe : evidence){
        if(pe.indelCount == 0){
            retainedEvidence.push_back(&pe);
        }
    }

    std::fprintf(stderr, "indel filter: %zu of %zu retained\n", retainedEvidence.size(), evidence.size());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(retainedEvidence.size());
        for(const amber::PositionEvidence *pe : retainedEvidence){
            rows.push_back({pe->chromosome, pe->position,
                    pe->chromosome + "\t" + std::to_string(pe->position)});
        }
        amber::CpDump::writeSorted("CP-A4", "chromosome\tposition", rows);
    }

    // ---- CP-A5：四道 filter 與排序 ----
    //
    // 對應 AmberApplication.runTumorOnly（AmberApplication.java:261-268）：
    //   ReadDepth >= TumorMinDepth（tumor-only 下為 25）
    //   aboveQualFilter（三個 filtered 計數器占 ReadDepth 的比例 < 0.15）
    //   RefSupport >= TumorOnlyMinSupport（2）
    //   AltSupport >= TumorOnlyMinSupport（2）
    // 之後 .sorted()——比較函式是 GenomePosition.compare，先比 ContigComparator 的
    // **rank 數值序**再比 position。rank 不是字串序：chr2 排在 chr10 之前。
    // idx 欄即此排序後的序位，驗收規則要求逐筆相同，故排序不可用字串比較。

    std::vector<const amber::PositionEvidence *> rawData;
    rawData.reserve(retainedEvidence.size());

    for(const amber::PositionEvidence *pe : retainedEvidence){
        if(pe->readDepth < amber::DEFAULT_TUMOR_ONLY_MIN_DEPTH){
            continue;
        }
        if(!amber::aboveQualFilter(*pe)){
            continue;
        }
        if(pe->refSupport < amber::DEFAULT_TUMOR_ONLY_MIN_SUPPORT){
            continue;
        }
        if(pe->altSupport < amber::DEFAULT_TUMOR_ONLY_MIN_SUPPORT){
            continue;
        }
        rawData.push_back(pe);
    }

    // Java 的 Stream.sorted() 是穩定排序；來源順序為 ArrayListMultimap 的走訪順序。
    // 位點的 (chromosome, position) 在此唯一，故排序為全序，穩定性不影響結果。
    std::stable_sort(rawData.begin(), rawData.end(),
            [](const amber::PositionEvidence *a, const amber::PositionEvidence *b){
                return amber::genomePositionLess(*a, *b);
            });

    std::fprintf(stderr, "four filters: %zu of %zu retained\n", rawData.size(), retainedEvidence.size());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(rawData.size());
        for(std::size_t i = 0; i < rawData.size(); ++i){
            const amber::PositionEvidence &pe = *rawData[i];
            std::string line = pe.chromosome;
            line += "\t" + std::to_string(pe.position);
            line += "\t" + std::to_string(i);
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
        // 排序本身即被比對，故用保留呼叫端順序的 write，不是 writeSorted
        amber::CpDump::write("CP-A5",
                "chromosome\tposition\tidx\tref\talt\treadDepth\tindelCount"
                "\trefSupport\taltSupport\tbaseQualFiltered\tmapQualFiltered\tseqTechFiltered",
                rows);
    }

    // ---- CP-A6：noise floor 與 contamination ----
    //
    // 對應 AmberApplication.runTumorOnly（AmberApplication.java:277-281）：
    // TumorOnlyPurityAnalysis 以 rawData（= CP-A5 的序列）為輸入，
    // cutoff() 取 min(MIN_CUTOFF, 最小的 copy-number peak / 3)，
    // contamination 取汙染 peak 的最大 vaf（無則 0.0）。

    std::vector<amber::PositionEvidence> rawDataValues;
    rawDataValues.reserve(rawData.size());
    for(const amber::PositionEvidence *pe : rawData){
        rawDataValues.push_back(*pe);
    }

    const amber::NoiseFloorResult noiseFloorResult = amber::computeNoiseFloor(rawDataValues);

    std::fprintf(stderr, "noise floor: %zu evidence points, %zu after immune filter, %zu maxima, "
            "noiseFloor(%.3f) contamination(%.3f)\n",
            noiseFloorResult.evidencePoints, noiseFloorResult.evidencePointsAfterImmuneFilter,
            noiseFloorResult.maximaDiagnostics.size(),
            noiseFloorResult.noiseFloor, noiseFloorResult.contamination);

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.push_back({"", 0, "noiseFloor\t" + amber::CpDump::num(noiseFloorResult.noiseFloor)});
        rows.push_back({"", 0, "contamination\t" + amber::CpDump::num(noiseFloorResult.contamination)});
        for(std::size_t i = 0; i < noiseFloorResult.contaminationPeakVafs.size(); ++i){
            rows.push_back({"", 0, "contaminationPeak." + std::to_string(i) + "\t"
                    + amber::CpDump::num(noiseFloorResult.contaminationPeakVafs[i])});
        }
        amber::CpDump::write("CP-A6", "field\tvalue", rows);

        // 診斷輸出（不在驗收規則內）：每個區域極大值的中間量，
        // 對應 Java 以 -log_debug 印出的同一組數字。CP-A6 只有三列，
        // 中間過程若有偏離不會在該檢查點顯現，故另存這一份供日後二分。
        std::vector<amber::CpDump::Row> diag;
        diag.push_back({"", 0, "evidencePoints\t" + std::to_string(noiseFloorResult.evidencePoints)});
        diag.push_back({"", 0, "afterImmuneFilter\t"
                + std::to_string(noiseFloorResult.evidencePointsAfterImmuneFilter)});
        for(const amber::PeakDiagnostic &d : noiseFloorResult.maximaDiagnostics){
            char buffer[256];
            std::snprintf(buffer, sizeof(buffer),
                    "peak.%.3f\tscore=%.6f;homProportion=%.6f;chrArmAUC=%.6f;mutationAUC=%.6f;captured=%d;%s",
                    d.vaf, d.score, d.homozygousProportion, d.chrArmAuc, d.mutationAuc,
                    d.capturedPoints, d.classification.c_str());
            diag.push_back({"", 0, buffer});
        }
        amber::CpDump::write("CP-A6-diagnostic", "field\tvalue", diag);
    }

    // ---- CP-A6b：noise floor 的套用 ----
    //
    // 對應 AmberApplication.java:284-288：兩個頻率都必須是有限值，
    // 且以 **Doubles.greaterOrEqual（epsilon 1e-10）** 而非 >= 與 noiseFloor 比較。
    // refFrequency / altFrequency 的分母是 **ReadDepth**，不是 Ref+Alt（TumorBAF.java:29-32）。

    std::vector<const amber::PositionEvidence *> tumorBAFList;
    tumorBAFList.reserve(rawData.size());

    for(const amber::PositionEvidence *pe : rawData){
        const double refFrequency = pe->refSupport / static_cast<double>(pe->readDepth);
        const double altFrequency = pe->altSupport / static_cast<double>(pe->readDepth);

        if(!std::isfinite(refFrequency)
                || !amber::doublesGreaterOrEqual(refFrequency, noiseFloorResult.noiseFloor)){
            continue;
        }
        if(!std::isfinite(altFrequency)
                || !amber::doublesGreaterOrEqual(altFrequency, noiseFloorResult.noiseFloor)){
            continue;
        }
        tumorBAFList.push_back(pe);
    }

    // Java 端在此再做一次 .sorted()；輸入已依同一比較函式排序且鍵唯一，故順序不變。
    std::stable_sort(tumorBAFList.begin(), tumorBAFList.end(),
            [](const amber::PositionEvidence *a, const amber::PositionEvidence *b){
                return amber::genomePositionLess(*a, *b);
            });

    std::fprintf(stderr, "noise floor applied: %zu of %zu retained\n", tumorBAFList.size(), rawData.size());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(tumorBAFList.size());
        for(std::size_t i = 0; i < tumorBAFList.size(); ++i){
            const amber::PositionEvidence &pe = *tumorBAFList[i];
            const double refFrequency = pe.refSupport / static_cast<double>(pe.readDepth);
            const double altFrequency = pe.altSupport / static_cast<double>(pe.readDepth);
            std::string line = pe.chromosome;
            line += "\t" + std::to_string(pe.position);
            line += "\t" + std::to_string(i);
            line += "\t" + amber::CpDump::num(refFrequency);
            line += "\t" + amber::CpDump::num(altFrequency);
            rows.push_back({pe.chromosome, pe.position, std::move(line)});
        }
        amber::CpDump::write("CP-A6b", "chromosome\tposition\tidx\trefFrequency\taltFrequency", rows);
    }

    // ---- CP-A7：AmberBAF 轉換 ----
    //
    // 對應 AmberApplication.java:290-291 與 AmberUtils.fromTumorBaf（AmberUtils.java:58-66）。
    // **tumorBAF 的分母是 Alt+Ref，不是 ReadDepth**；tumorDepth 才是 ReadDepth。
    // tumor-only 下 NormalAltSupport / NormalRefSupport / NormalReadDepth 皆為 0
    // （TumorBAF.fromNormal 從歸零的 PositionEvidence 複製），故 normalBAF = 0/0 = NaN。

    std::vector<amber::AmberBAF> amberBAFList;
    amberBAFList.reserve(tumorBAFList.size());

    for(const amber::PositionEvidence *pe : tumorBAFList){
        const int tumorAltCount = pe->altSupport;
        const double tumorBaf = tumorAltCount / static_cast<double>(tumorAltCount + pe->refSupport);

        const int normalAltCount = 0;
        const int normalRefSupport = 0;
        const double normalBaf = normalAltCount / static_cast<double>(normalAltCount + normalRefSupport);

        amber::AmberBAF baf;
        baf.chromosome = pe->chromosome;
        baf.position = pe->position;
        baf.tumorBAF = tumorBaf;
        baf.tumorDepth = pe->readDepth;
        baf.normalBAF = normalBaf;
        baf.normalDepth = 0;

        // AmberApplication.java:291 的 filter(x -> Double.isFinite(x.tumorBAF()))
        if(!std::isfinite(baf.tumorBAF)){
            continue;
        }

        amberBAFList.push_back(std::move(baf));
    }

    std::fprintf(stderr, "amber BAF: %zu of %zu retained\n", amberBAFList.size(), tumorBAFList.size());

    if(amber::CpDump::enabled()){
        std::vector<amber::CpDump::Row> rows;
        rows.reserve(amberBAFList.size());
        for(std::size_t i = 0; i < amberBAFList.size(); ++i){
            const amber::AmberBAF &b = amberBAFList[i];
            std::string line = b.chromosome;
            line += "\t" + std::to_string(b.position);
            line += "\t" + std::to_string(i);
            line += "\t" + amber::CpDump::num(b.tumorBAF);
            line += "\t" + amber::CpDump::num(b.tumorModifiedBAF());
            line += "\t" + std::to_string(b.tumorDepth);
            line += "\t" + amber::CpDump::num(b.normalBAF);
            line += "\t" + std::to_string(b.normalDepth);
            rows.push_back({b.chromosome, b.position, std::move(line)});
        }
        amber::CpDump::write("CP-A7",
                "chromosome\tposition\tidx\ttumorBAF\ttumorModifiedBAF"
                "\ttumorDepth\tnormalBAF\tnormalDepth", rows);
    }

    // ---- stage 輸出：amber.baf.tsv.gz 與 amber.qc ----
    //
    // 對應 ResultsWriter.persistBAF / persistQC（ResultsWriter.java:38-61）。
    // persistBAF 內另有 PCF 分段（同檔 43-51），屬 EXP-010 的範圍，此處不實作。
    if(!outputDir.empty()){
        if(sampleId.empty()){
            throw std::runtime_error("-output_dir 需要同時給定 -tumor");
        }

        const std::string bafPath = outputDir + "/" + sampleId + ".amber.baf.tsv.gz";
        amber::writeAmberBafFile(bafPath, amberBAFList);

        // AmberApplication.java:293 persistQC(0, contamination, null)：
        // tumor-only 下 consanguinityProportion 為 0、uniparentalDisomy 為 null
        const std::string qcPath = outputDir + "/" + sampleId + ".amber.qc";
        amber::writeAmberQcFile(qcPath, noiseFloorResult.contamination, 0.0);

        std::fprintf(stderr, "wrote %s and %s\n", bafPath.c_str(), qcPath.c_str());
    }

    // ---- CP-A8 / CP-A9：PCF 分段 ----
    //
    // 對應 ResultsWriter.persistBAF 內的分段（ResultsWriter.java:43-51）→
    // BAFSegmenter.writeSegments → PerArmSegmenter。gamma 硬編碼 100.0，AMBER 4.3 無 CLI 可調。
    // 注意 tumor-only 也會做分段：runTumorOnly 本身沒呼叫，但它呼叫的 persistBAF 內有。

    const amber::SegmentationResult segmentation = amber::segmentBafs(amberBAFList);

    std::size_t segmentCount = 0;
    for(const amber::ArmSegments &arm : segmentation.arms){
        segmentCount += arm.segments.size();
    }
    std::fprintf(stderr, "PCF segmentation: %d values across %zu arms, penaltyMode(%s), %zu segments\n",
            segmentation.totalCount, segmentation.arms.size(),
            segmentation.penaltyMode.c_str(), segmentCount);

    amber::writeSegmentationCheckpoints(segmentation);

    if(!outputDir.empty() && !sampleId.empty()){
        const std::string pcfPath = outputDir + "/" + sampleId + ".amber.baf.pcf";
        amber::writeSegmentsFile(pcfPath, segmentation);
        std::fprintf(stderr, "wrote %s\n", pcfPath.c_str());
    }

    return 0;
}
