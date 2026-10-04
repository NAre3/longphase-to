#ifndef AMBER_POSITIONEVIDENCE_H
#define AMBER_POSITIONEVIDENCE_H

#include <string>

namespace amber {

// 對應 amber-v4.3 的 PositionEvidence（PositionEvidence.java:12-41）。
// 七個計數器的語義見 behaviour-contract.md §2.1。
struct PositionEvidence
{
    std::string chromosome;
    int position = 0;
    char ref = 'N';
    char alt = 'N';

    int readDepth = 0;
    int indelCount = 0;
    int refSupport = 0;
    int altSupport = 0;
    int baseQualFiltered = 0;
    int mapQualFiltered = 0;
    int seqTechFiltered = 0;
};

}

#endif
