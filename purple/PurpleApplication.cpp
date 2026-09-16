#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include "PurpleInput.h"
#include "PurpleSegmentation.h"
#include "PurpleObserved.h"
#include "PurpleFitting.h"
#include "PurpleCopyNumber.h"
#include "PurpleSummary.h"
#include "PurpleWriters.h"
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
        const std::string ensembl = value(argc, argv, "-ensembl_data_dir");
        const std::string output = value(argc, argv, "-output_dir");
        if(sample.empty() || amber.empty() || cobalt.empty() || reference.empty() || ensembl.empty()){
            std::cerr << "usage: purple_port -tumor <sample> -amber <dir> -cobalt <dir> -ref_genome <fasta> -ensembl_data_dir <dir> [-output_dir <dir>] [-cpdump_dir <dir>]\n";
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
        const auto bestFit = purple::selectTumorOnlyBestFit(fits, observed);
        purple::dumpBestFit(bestFit);
        const auto fittedRegions = purple::fitObservedRegions(observed, bestFit.fit, inputs.averageTumorDepth, inputs.cobaltGender);
        purple::dumpFittedRegions(fittedRegions);
        const auto copyNumbers = purple::buildCopyNumbers(fittedRegions, bestFit.fit, inputs.cobaltGender);
        purple::dumpCopyNumbers(copyNumbers);
        const auto summary = purple::buildSummaryContext(inputs, bestFit, copyNumbers, ensembl);
        purple::dumpSummaryContext(inputs, bestFit, copyNumbers, summary);
        purple::writeCoreOutputs(output, inputs, fits, bestFit, fittedRegions, copyNumbers, summary);
        std::cerr << "PURPLE P10 core-output stage complete: " << inputs.bafs.size() << " BAF, "
                  << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
                  << inputs.cobaltTumorPcf.size() << " Cobalt PCF, " << segments.size() << " support segments, "
                  << observed.size() << " observed regions, " << fits.size() << " purity/ploidy candidates\n";
        return 0;
    }catch(const std::exception &exception){
        std::cerr << "purple_port: " << exception.what() << '\n';
        return 1;
    }
}
