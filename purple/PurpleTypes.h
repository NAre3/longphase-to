#ifndef PURPLE_TYPES_H
#define PURPLE_TYPES_H

#include <cstdint>
#include <string>
#include <vector>

namespace purple {

enum class Gender { FEMALE, MALE };

inline const char *genderName(Gender gender){ return gender == Gender::FEMALE ? "FEMALE" : "MALE"; }

struct AmberBaf {
    std::string chromosome;
    int position = 0;
    double tumorBaf = 0;
    int tumorDepth = 0;
    double normalBaf = 0;
    int normalDepth = 0;
};

struct CobaltRatio {
    std::string chromosome;
    int position = 0;
    double referenceReadDepth = 0;
    double referenceGcRatio = 0;
    double referenceGcContent = 0;
    double referenceGcDiploidRatio = 0;
    double tumorReadDepth = 0;
    double tumorGcRatio = 0;
    double tumorGcContent = 0;
};

enum class PcfSource { REFERENCE_RATIO, TUMOR_RATIO, TUMOR_BAF };

inline const char *pcfSourceName(PcfSource source){
    if(source == PcfSource::REFERENCE_RATIO){ return "REFERENCE_RATIO"; }
    if(source == PcfSource::TUMOR_RATIO){ return "TUMOR_RATIO"; }
    return "TUMOR_BAF";
}

struct PcfPosition {
    PcfSource source = PcfSource::TUMOR_RATIO;
    std::string chromosome;
    int position = 0;
    int minPosition = 0;
    int maxPosition = 0;
    bool isSegmentEnd() const { return position == minPosition; }
    bool isSegmentStart() const { return !isSegmentEnd(); }
};

struct InputData {
    std::string sampleId;
    std::vector<AmberBaf> bafs;
    std::vector<CobaltRatio> ratios;
    std::vector<PcfPosition> amberPcf;
    std::vector<PcfPosition> cobaltTumorPcf;
    double contamination = 0;
    int averageTumorDepth = 100;
    Gender amberGender = Gender::MALE;
    Gender cobaltGender = Gender::MALE;
};

enum class SegmentSupport { NONE, TELOMERE, CENTROMERE, MULTIPLE, EXCL };

inline const char *segmentSupportName(SegmentSupport support){
    if(support == SegmentSupport::TELOMERE){ return "TELOMERE"; }
    if(support == SegmentSupport::CENTROMERE){ return "CENTROMERE"; }
    if(support == SegmentSupport::MULTIPLE){ return "MULTIPLE"; }
    if(support == SegmentSupport::EXCL){ return "EXCL"; }
    return "NONE";
}

struct SupportSegment {
    std::string chromosome;
    int start = 0;
    int end = 0;
    bool ratioSupport = true;
    SegmentSupport support = SegmentSupport::NONE;
    bool svCluster = false;
    int minStart = 0;
    int maxStart = 0;
};

enum class GermlineStatus { UNKNOWN, DIPLOID, CENTROMETIC, EXCLUDED };

inline const char *germlineStatusName(GermlineStatus status){
    if(status == GermlineStatus::DIPLOID){ return "DIPLOID"; }
    if(status == GermlineStatus::CENTROMETIC){ return "CENTROMETIC"; }
    if(status == GermlineStatus::EXCLUDED){ return "EXCLUDED"; }
    return "UNKNOWN";
}

struct ObservedRegion {
    SupportSegment segment;
    int bafCount = 0;
    double observedBaf = 0;
    int depthWindowCount = 0;
    double observedTumorRatio = 0;
    double observedNormalRatio = 0;
    double unnormalisedObservedNormalRatio = 0;
    GermlineStatus germlineStatus = GermlineStatus::UNKNOWN;
    double gcContent = 0;
    double minorAlleleCopyNumberDeviation = 0;
    double majorAlleleCopyNumberDeviation = 0;
    double deviationPenalty = 0;
    double eventPenalty = 0;
    double refNormalisedCopyNumber = 0;
    double tumorCopyNumber = 0;
    double tumorBaf = 0;
    double fittedTumorCopyNumber = 0;
    double fittedBaf = 0;
};

struct FittedPurity {
    double purity = 0;
    double normFactor = 0;
    double ploidy = 0;
    double score = 0;
    double diploidProportion = 0;
    double somaticPenalty = 0;
};

struct FittedPurityScore {
    double minPurity = 0, maxPurity = 0;
    double minPloidy = 0, maxPloidy = 0;
    double minDiploidProportion = 0, maxDiploidProportion = 0;
};

struct BestFit {
    FittedPurity fit;
    FittedPurityScore score;
    std::string method;
};

enum class CopyNumberMethod { UNKNOWN, BAF_WEIGHTED, LONG_ARM };
inline const char *copyNumberMethodName(CopyNumberMethod method){
    if(method == CopyNumberMethod::BAF_WEIGHTED){ return "BAF_WEIGHTED"; }
    if(method == CopyNumberMethod::LONG_ARM){ return "LONG_ARM"; }
    return "UNKNOWN";
}

struct PurpleCopyNumber {
    std::string chromosome;
    int start=0,end=0,bafCount=0,depthWindowCount=0,minStart=0,maxStart=0;
    double averageActualBaf=0,averageObservedBaf=0,averageTumorCopyNumber=0,gcContent=0;
    SegmentSupport segmentStartSupport=SegmentSupport::NONE,segmentEndSupport=SegmentSupport::NONE;
    CopyNumberMethod method=CopyNumberMethod::UNKNOWN;
};

}

#endif
