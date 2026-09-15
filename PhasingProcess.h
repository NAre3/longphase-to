#ifndef PHASINGPROCESS_H
#define PHASINGPROCESS_H
#include "MethylXgbModel.h"
#include "Util.h"

struct PhasingParameters
{
    int numThreads;
    int distance;
    std::string snpFile;
    std::string svFile;
    std::vector<std::string> bamFile;
    std::string modFile="";
    std::string ponFile="";
    std::string strictPonFile="";
    std::string fastaFile;
    std::string resultPrefix;
    std::string callerStr;
    Caller caller;
    bool generateDot;
    bool phaseIndel;
    bool disablePonTag=true;
    bool disableCalling=false;
    bool disableRefineSomatic=false;
    bool outputLOH = false;  // Whether to output LOH results
    bool outputSGE = false;  // Whether to output SmallGenomicEvent results
    bool outputLGE = false;  // Whether to output LargeGenomicEvent results
    bool outputGE = false;   // Whether to output GenomicEvent results
    bool enableMethylXgb = true;
    
    int connectAdjacent;
    int mappingQuality;
    double mismatchRate;
    
    int baseQuality;
    double edgeWeight;
    
    double snpConfidence;
    double readConfidence;
    
    double edgeThreshold;
    double overlapThreshold;
    
    int somaticConnectAdjacent;
    int methylXgbWindow = 2000;
    float methylXgbMethHigh = 0.8f;
    float methylXgbMethLow = 0.2f;
    double methylXgbSnvThreshold = METHYL_XGB_DEFAULT_SNV_THRESHOLD;
    double methylXgbIndelThreshold = METHYL_XGB_DEFAULT_INDEL_THRESHOLD;
    
    std::string version;
    std::string command;
    
    // If negative, purity is not provided by user and should be estimated
    double purity = -1.0;

    // ---- AMBER/COBALT 整合（EXP-I02，D-I2）----
    // 四個旗標皆選填。**全部未給定 → 完全不啟用**，LongPhase-TO 的行為與整合前一致
    // （這是 F3 迴歸防線的對照基準）。給了一部分而缺另一部分即報錯，不默默半啟用。
    // 本票（EXP-I02）只接 AMBER 的兩個；gcProfile / diploidBed 於 EXP-I03 接上。
    std::string amberLoci = "";
    std::string amberExcludedBed = "";
    std::string amberOutputDir = "";
    std::string amberSampleId = "";
    int amberMinBaseQuality = 13;   // amber::DEFAULT_MIN_BASE_QUALITY
    int amberMinMapQuality = 50;    // amber::DEFAULT_MIN_MAPPING_QUALITY
    std::string amberCpDumpDir = "";

    bool amberEnabled() const { return !amberLoci.empty(); }
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif
