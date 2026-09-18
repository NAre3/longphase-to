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

    // EXP-I03：COBALT 的三個必填資源 + 選用輸出。
    std::string cobaltGcProfile = "";
    std::string cobaltDiploidBed = "";
    std::string cobaltExcludedRegions = "";
    std::string cobaltOutputDir = "";
    std::string cobaltSampleId = "";
    int cobaltMinMapQuality = 10;   // cobalt::DEFAULT -min_quality

    // PURPLE 整合。PURPLE **不讀 BAM**，消費的是 AMBER/COBALT 的結果，
    // 因此它接在兩者的 postscan 之後，不參與共用掃描。
    //
    // 刻意**沒有** ref_genome 欄位：purple_port 用 -ref_genome 只為了讀 .fai 取
    // 染色體長度，而此處長度直接取自 cobaltPrescan.chromosomes（來源是 BAM
    // header @SQ）。實測 HCC1937_t50_n00 的 header 195 條 contig 長度與 .fai
    // 逐一相同，故不必再要求使用者給一次參考基因體。
    //
    // 中間檔案預設不落地：AMBER/COBALT 的 stage 輸出改以記憶體交給 PURPLE，
    // 只有在使用者明確給了 amberOutputDir / cobaltOutputDir 時才會寫出。
    std::string purpleEnsemblDataDir = "";
    std::string purpleOutputDir = "";
    std::string purpleSampleId = "";
    std::string purpleCpDumpDir = "";

    bool amberEnabled() const { return !amberLoci.empty(); }
    bool cobaltEnabled() const { return !cobaltGcProfile.empty(); }
    // PURPLE 需要 AMBER 與 COBALT 兩邊的結果，缺一不可。
    bool purpleEnabled() const {
        return !purpleEnsemblDataDir.empty() && amberEnabled() && cobaltEnabled();
    }
};

class PhasingProcess
{

    public:
        PhasingProcess(PhasingParameters params);
        ~PhasingProcess();

};


#endif
