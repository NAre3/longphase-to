// AMBER 4.3 tumor-only 移植：amber_port 的進入點。
//
// EXP-I02（共用掃描層整合）把原本寫在此檔 main() 內的三段流程搬到 AmberPipeline.{h,cpp}：
//   prescan（CP-A1→CP-A2b）／BAM 掃描／postscan（CP-A3→CP-A9 與三個 stage 輸出）。
// 本檔自此只剩「解析 CLI」與「照順序呼叫那三段」，**沒有任何計算**。
//
// 整合進 LongPhase-TO 的版本呼叫的是同一組 prescan/postscan，差別只在中間那段掃描
// 由共用走訪的 amber::ContigSink 取代 amber::processBam。
// 這是「凍結候選（91959e4 / b71e68e）的行為未因整合而改變」這句話的依據。
//
// 逐一對齊的行為出處：RUN-004/behaviour-contract.md

#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "AmberPipeline.h"
#include "BamEvidenceReader.h"
#include "../common/CpDump.h"

namespace {

// 兩個保真度驗證用的旗標**刻意不列進 usage**（owner 裁示 2026-09-02）：
//
//   -cpdump_dir <dir>     開啟 checkpoint dump。未給定時 CpDump::enabled() 為 false，
//                         所有 dump 呼叫直接跳過，不產生檔案、不影響計算。
//   -debug_only_chr <chr> 只保留該染色體的 loci，用於縮短迭代週期。
//                         **不是移植的行為**——AMBER 的 -specific_chr 並不限制 loci。
//                         正式記錄的執行必須不帶此旗標。
//
// 兩者都是「怎麼驗證這份移植」的基礎設施，不是 AMBER 的功能，因此不對一般使用者呈現；
// 但保留在程式碼裡，因為重驗的需求會實際發生（例如換用帶 Frequency 欄的 loci 檔時，
// 見 AmberSitesFile.cpp 的 F-R1 註解）。用法見
// research/studies/purple-port-amber-fidelity-v1/ 的各 run 執行腳本。
std::string usage(){
    return "usage: amber_port -loci <AmberGermlineSites.tsv.gz> "
           "-tumor_only_excluded_bed <tumorOnlyExcludedSnp.38.bed>\n"
           "                  [-tumor_bam <bam>] [-min_base_quality N] [-min_map_quality N]\n"
           "                  [-output_dir <dir>] [-tumor <sampleId>] [-threads N]\n"
           "                  [-write_tumor_data] [-write_version]\n"
           "\n"
           "-tumor_bam 未給定時只跑到 CP-A2 為止。\n"
           "-output_dir 與 -tumor 同時給定時寫出 amber.baf.tsv.gz、amber.qc 與 amber.baf.pcf\n"
           "  ——這三個即 PURPLE 唯一會讀取的檔案（purple/AmberData.java）。\n"
           "-write_tumor_data 另外寫出 <sample>.amber.tumor.raw.tsv.gz（預設關閉）。\n"
           "-write_version 另外寫出 amber.version（預設關閉）。\n"
           "  兩者皆為稽核／除錯用，PURPLE 不讀取，且已實證不影響任何計算。\n";
}

}

int main(int argc, char **argv){
    amber::PipelineConfig cfg;
    std::string cpDumpDir;

    // 以下三個區域變數是為了讓搬移前的參數解析程式碼逐字沿用，解析完再填回 cfg。
    std::string lociPath;
    std::string bedPath;
    std::string tumorBam;
    std::string debugOnlyChr;
    std::string outputDir;
    int threads = 1;
    bool writeTumorData = false;
    bool writeVersion = false;
    std::string sampleId;
    int minBaseQuality = amber::DEFAULT_MIN_BASE_QUALITY;
    int minMappingQuality = amber::DEFAULT_MIN_MAPPING_QUALITY;

    // 注意：迴圈上界必須是 argc 而非 argc-1。原本寫成 argc-1 是為「旗標＋值」成對設計，
    // 但那會讓**位於最後一個位置的布林旗標**被靜默忽略（-write_version 曾因此無效）。
    // 取值型參數改為明確檢查後面還有沒有東西。
    const auto needValue = [&](int i, const char *name) -> const char * {
        if(i + 1 >= argc){
            throw std::runtime_error(std::string(name) + " 需要一個值");
        }
        return argv[i + 1];
    };

    for(int i = 1; i < argc; ++i){
        const std::string arg = argv[i];
        if(arg == "-loci"){
            lociPath = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-tumor_only_excluded_bed"){
            bedPath = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-cpdump_dir"){
            cpDumpDir = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-tumor_bam"){
            tumorBam = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-min_base_quality"){
            minBaseQuality = std::stoi(needValue(i, arg.c_str())); ++i;
        }else if(arg == "-min_map_quality"){
            minMappingQuality = std::stoi(needValue(i, arg.c_str())); ++i;
        }else if(arg == "-debug_only_chr"){
            debugOnlyChr = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-write_tumor_data"){
            writeTumorData = true;
        }else if(arg == "-write_version"){
            writeVersion = true;
        }else if(arg == "-threads"){
            threads = std::stoi(needValue(i, arg.c_str())); ++i;
        }else if(arg == "-output_dir"){
            outputDir = needValue(i, arg.c_str()); ++i;
        }else if(arg == "-tumor"){
            sampleId = needValue(i, arg.c_str()); ++i;
        }
    }

    if(lociPath.empty() || bedPath.empty()){
        std::cerr << usage();
        return 1;
    }

    cfg.lociPath = lociPath;
    cfg.bedPath = bedPath;
    cfg.tumorBam = tumorBam;
    cfg.debugOnlyChr = debugOnlyChr;
    cfg.outputDir = outputDir;
    cfg.sampleId = sampleId;
    cfg.threads = threads;
    cfg.writeTumorData = writeTumorData;
    cfg.writeVersion = writeVersion;
    cfg.minBaseQuality = minBaseQuality;
    cfg.minMappingQuality = minMappingQuality;

    lp::CpDump::setDir(cpDumpDir);

    amber::PrescanResult prescan = amber::prescan(cfg);

    // tumorBam 未給定時只跑到 CP-A2 為止（usage 明列的行為）
    if(!prescan.hasBamStage){
        return 0;
    }

    const amber::BamScanStats stats = amber::processBam(
            cfg.tumorBam, prescan.tasks, cfg.minMappingQuality, cfg.minBaseQuality, cfg.threads);

    amber::postscan(cfg, prescan.evidence, stats);

    return 0;
}
