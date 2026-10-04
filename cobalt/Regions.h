#ifndef COBALT_REGIONS_H
#define COBALT_REGIONS_H

#include <string>
#include <vector>

namespace cobalt {

// ExcludedRegionsFile：jar 內部資源 /regions/excluded_regions_v38.tsv（G8）。
// 染色體名無 chr 前綴；順序 = 檔案順序。
struct ExcludedRegion
{
    std::string chromosome;
    int start = 0;
    int end = 0;
};

std::vector<ExcludedRegion> loadExcludedRegions(const std::string &path);

// DiploidRegionLoader 產生的逐 window 狀態。
struct DiploidStatus
{
    std::string contig;
    int start = 0;
    int end = 0;
    bool isDiploid = false;
};

struct DiploidRegions
{
    std::vector<std::string> chromosomes;                 // 短名，依 ordinal 昇冪
    std::vector<std::vector<DiploidStatus>> byChromosome;

    bool empty() const { return chromosomes.empty(); }
    std::size_t totalEntries() const;
    const std::vector<DiploidStatus> *find(const std::string &shortName) const;
};

DiploidRegions loadDiploidRegions(const std::string &bedPath);

}

#endif
