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
            std::cerr << "usage: purple_port -tumor <sample> -amber <dir> -cobalt <dir> -ref_genome <fasta> -ensembl_data_dir <dir> [-output_dir <dir>] [-cpdump_dir <dir>]\n";
            return 2;
        }
        cfg.threads = std::max(1, std::atoi(value(argc, argv, "-threads").c_str()));
        lp::CpDump::setDir(value(argc, argv, "-cpdump_dir"));
        purple::run(cfg);
        return 0;
    }catch(const std::exception &exception){
        std::cerr << "purple_port: " << exception.what() << '\n';
        return 1;
    }
}
