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

}

#endif
