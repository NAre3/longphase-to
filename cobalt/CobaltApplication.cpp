// COBALT 移植：cobalt_port 的進入點。
//
// EXP-I03（共用掃描層整合）把原本寫在此檔 main() 內的三段流程搬到 CobaltPipeline.{h,cpp}：
//   prescan（CP-C1→CP-C5）／BAM 掃描／postscan（CP-C6→CP-C14 與三個 stage 輸出）。
// 本檔自此只剩「解析 CLI」與「照順序呼叫那三段」，**沒有任何計算**。
//
// 整合進 LongPhase-TO 的版本呼叫的是同一組 prescan/postscan，差別只在中間那段掃描
// 由共用走訪的 cobalt::DepthSink 取代 cobalt::calculateReadDepths。
// 這是「凍結候選（b71e68e）的行為未因整合而改變」這句話的依據。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "CobaltPipeline.h"
#include "ReadDepth.h"
#include "../common/CpDump.h"

using lp::CpDump;

namespace {

std::string arg(int argc, char **argv, const char *flag, const std::string &fallback = std::string())
{
    for(int i = 1; i + 1 < argc; ++i)
    {
        if(std::strcmp(argv[i], flag) == 0){ return argv[i + 1]; }
    }
    return fallback;
}

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

    cobalt::PipelineConfig cfg;
    cfg.bamPath = bamPath;
    cfg.gcProfile = gcProfile;
    cfg.diploidBed = diploidBed;
    cfg.excludedPath = excludedPath;
    cfg.outputDir = arg(argc, argv, "-output_dir", ".");
    cfg.sampleId = arg(argc, argv, "-tumor", "tumor");
    cfg.threads = threads;
    cfg.minMappingQuality = minMappingQuality;
    cfg.pcfGamma = std::atof(arg(argc, argv, "-pcf_gamma", "100").c_str());
    cfg.includeDuplicates = includeDuplicates;

    const cobalt::PrescanResult pre = cobalt::prescan(cfg);

    const std::vector<cobalt::DepthReading> depths = cobalt::calculateReadDepths(
            cfg.bamPath, pre.chromosomes, pre.partitions,
            cfg.minMappingQuality, cfg.includeDuplicates, cfg.threads);

    cobalt::postscan(cfg, pre, depths);

    return 0;
}
