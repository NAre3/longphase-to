#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

#include "PurpleInput.h"
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
        if(sample.empty() || amber.empty() || cobalt.empty()){
            std::cerr << "usage: purple_port -tumor <sample> -amber <dir> -cobalt <dir> [-cpdump_dir <dir>]\n";
            return 2;
        }
        lp::CpDump::setDir(value(argc, argv, "-cpdump_dir"));
        const purple::InputData inputs = purple::loadTumorOnlyInputs(sample, amber, cobalt);
        purple::dumpInputCheckpoint(inputs);
        std::cerr << "PURPLE P02 input stage complete: " << inputs.bafs.size() << " BAF, "
                  << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
                  << inputs.cobaltTumorPcf.size() << " Cobalt PCF\n";
        return 0;
    }catch(const std::exception &exception){
        std::cerr << "purple_port: " << exception.what() << '\n';
        return 1;
    }
}
