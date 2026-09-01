// AMBER 4.3 tumor-only 移植：第一段（CP-A1 site 載入、CP-A2 tumor-only blacklist 過濾）
//
// 對照對象是 hmftools tag amber-v4.3 的 AmberApplication.loadAmberSites()
// 與 AmberApplication.hetLociTumorOnly()。本階段刻意不共用 LongPhase-TO 既有的
// BAM 掃描，也不接進主程式的建置，理由見研究規格 §0。

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "AmberSitesFile.h"
#include "ChrBaseRegion.h"
#include "CpDump.h"
#include "HumanChromosome.h"

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
           "-tumor_only_excluded_bed <tumorOnlyExcludedSnp.38.bed> -cpdump_dir <dir>\n";
}

}

int main(int argc, char **argv){
    std::string lociPath;
    std::string bedPath;
    std::string cpDumpDir;

    for(int i = 1; i < argc - 1; ++i){
        const std::string arg = argv[i];
        if(arg == "-loci"){
            lociPath = argv[++i];
        }else if(arg == "-tumor_only_excluded_bed"){
            bedPath = argv[++i];
        }else if(arg == "-cpdump_dir"){
            cpDumpDir = argv[++i];
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

    return 0;
}
