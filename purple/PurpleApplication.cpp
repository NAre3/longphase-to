#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include "PurpleInput.h"
#include "PurpleSegmentation.h"
#include "PurpleObserved.h"
#include "PurpleFitting.h"
#include "../common/CpDump.h"

namespace {

std::string value(int argc, char **argv, const char *flag){
    for(int i = 1; i + 1 < argc; ++i){ if(std::strcmp(argv[i], flag) == 0){ return argv[i + 1]; } }
    return {};
}

}

int main(int argc, char **argv){
    try{
        const std::string sample = value(argc, argv, "-tumor");
        const std::string amber = value(argc, argv, "-amber");
        const std::string cobalt = value(argc, argv, "-cobalt");
        const std::string reference = value(argc, argv, "-ref_genome");
        if(sample.empty() || amber.empty() || cobalt.empty() || reference.empty()){
            std::cerr << "usage: purple_port -tumor <sample> -amber <dir> -cobalt <dir> -ref_genome <fasta> [-cpdump_dir <dir>]\n";
            return 2;
        }
        lp::CpDump::setDir(value(argc, argv, "-cpdump_dir"));
        const purple::InputData inputs = purple::loadTumorOnlyInputs(sample, amber, cobalt);
        purple::dumpInputCheckpoint(inputs);
        const auto segments = purple::createSupportSegments(inputs, reference);
        purple::dumpSupportSegments(segments);
        const auto observed = purple::createObservedRegions(inputs, segments);
        purple::dumpObservedRegions(observed);
        const auto fittingRegions = purple::selectFittingRegions(observed);
        purple::dumpFittingRegions(fittingRegions);
        const int threads = std::max(1, std::atoi(value(argc, argv, "-threads").c_str()));
        const auto fits = purple::fitPurityGrid(fittingRegions, inputs.averageTumorDepth, threads);
        purple::dumpPurityGrid(fits);
        std::cerr << "PURPLE P05 fitting-grid stage complete: " << inputs.bafs.size() << " BAF, "
                  << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
                  << inputs.cobaltTumorPcf.size() << " Cobalt PCF, " << segments.size() << " support segments, "
                  << observed.size() << " observed regions, " << fits.size() << " purity/ploidy candidates\n";
        return 0;
    }catch(const std::exception &exception){
        std::cerr << "purple_port: " << exception.what() << '\n';
        return 1;
    }
}
