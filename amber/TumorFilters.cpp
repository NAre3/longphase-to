#include "TumorFilters.h"

#include <cctype>
#include <cstdlib>

#include "HumanChromosome.h"

namespace amber {

namespace {

bool equalsIgnoreCase(const std::string &a, const char *b)
{
    std::size_t i = 0;
    for(; a[i] != '\0' && b[i] != '\0'; ++i){
        if(std::toupper(static_cast<unsigned char>(a[i])) != std::toupper(static_cast<unsigned char>(b[i]))){
            return false;
        }
    }
    return i == a.size() && b[i] == '\0';
}

}

bool aboveQualFilter(const PositionEvidence &pe)
{
    if(pe.readDepth == 0){
        return false;
    }

    const int filteredCount = pe.baseQualFiltered + pe.mapQualFiltered + pe.seqTechFiltered;
    const double filteredPercent = filteredCount / static_cast<double>(pe.readDepth);
    return filteredPercent < QUAL_FILTERED_THRESHOLD;
}

int chromosomeRank(const std::string &chromosome)
{
    const std::string trimmed = stripChrPrefix(chromosome);

    if(equalsIgnoreCase(trimmed, "X")){
        return 23;
    }
    if(equalsIgnoreCase(trimmed, "Y")){
        return 24;
    }
    if(equalsIgnoreCase(trimmed, "MT") || equalsIgnoreCase(trimmed, "M")){
        return 25;
    }

    // 對應 Integer.parseInt + NumberFormatException → INVALID_CHR_RANK
    if(trimmed.empty()){
        return 26;
    }

    std::size_t consumed = 0;
    long value = 0;
    try{
        value = std::stol(trimmed, &consumed);
    }catch(...){
        return 26;
    }

    if(consumed != trimmed.size()){
        return 26;
    }

    return static_cast<int>(value);
}

bool genomePositionLess(const PositionEvidence &a, const PositionEvidence &b)
{
    const int rankA = chromosomeRank(a.chromosome);
    const int rankB = chromosomeRank(b.chromosome);

    if(rankA != rankB){
        return rankA < rankB;
    }

    return a.position < b.position;
}

}
