#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include "PurplePipeline.h"
#include "../common/CpDump.h"

namespace {

std::string value(int argc, char **argv, const char *flag){
    for(int i = 1; i + 1 < argc; ++i){ if(std::strcmp(argv[i], flag) == 0){ return argv[i + 1]; } }
    return {};
}

}

// purple_port 的 main：引數解析 + PurplePipeline::run 的薄包裝。
// 流程主體在 PurplePipeline.cpp，與 longphase-to 整合版呼叫同一組函式。
// ［2026-10-03 註：整合版預設開啟 LOH 補回；purple_port 只在給 -loh_bed 時補，見 PurplePipeline.h。］
int main(int argc, char **argv){
    try{
        purple::PipelineConfig cfg;
        cfg.sampleId = value(argc, argv, "-tumor");
        cfg.amberDir = value(argc, argv, "-amber");
        cfg.cobaltDir = value(argc, argv, "-cobalt");
        cfg.refGenome = value(argc, argv, "-ref_genome");
        cfg.ensemblDataDir = value(argc, argv, "-ensembl_data_dir");
        cfg.outputDir = value(argc, argv, "-output_dir");
        if(cfg.sampleId.empty() || cfg.amberDir.empty() || cfg.cobaltDir.empty()
                || cfg.refGenome.empty() || cfg.ensemblDataDir.empty()){
            std::cerr << "usage: purple_port -tumor <sample> -amber <dir> -cobalt <dir> -ref_genome <fasta> -ensembl_data_dir <dir> [-output_dir <dir>] [-cpdump_dir <dir>] [-loh_bed <longphase-to _LOH.bed> -amber_loci <AmberGermlineSites.tsv.gz>]\n"
                      << "  -loh_bed and -amber_loci together enable the LOH fill (PurpleLohFill.h); without them purple_port does not fill.\n";
            return 2;
        }
        cfg.threads = std::max(1, std::atoi(value(argc, argv, "-threads").c_str()));
        lp::CpDump::setDir(value(argc, argv, "-cpdump_dir"));
        // 整合版的 LOH 來自記憶體；purple_port 沒有 LOH 來源，只在明確給了 BED 時才補。
        const std::string lohBed = value(argc, argv, "-loh_bed");
        const std::string amberLoci = value(argc, argv, "-amber_loci");
        if(lohBed.empty() != amberLoci.empty()){
            std::cerr << "purple_port: -loh_bed and -amber_loci must be given together\n";
            return 2;
        }
        if(!lohBed.empty()){
            purple::LohFillInput lohFill;
            lohFill.lohSpans = purple::loadLohBed(lohBed);
            lohFill.panelPositions = purple::loadPanelPositions(amberLoci);
            purple::run(cfg, &lohFill);
        }else{
            purple::run(cfg);
        }
        return 0;
    }catch(const std::exception &exception){
        std::cerr << "purple_port: " << exception.what() << '\n';
        return 1;
    }
}
