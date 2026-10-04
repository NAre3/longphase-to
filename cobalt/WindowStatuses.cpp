#include "WindowStatuses.h"

#include <algorithm>
#include <set>
#include <utility>

#include "CobaltConstants.h"
#include "../common/HumanChromosome.h"

namespace cobalt {

namespace {

// WindowStatuses.indexFor：整數除法，0-indexed（start = 1 → 0）
inline int indexFor(int position){ return position / WINDOW_SIZE; }

// ChrBaseRegion.overlaps：兩區間相交（皆 1-based 含端點）
inline bool overlaps(int aStart, int aEnd, int bStart, int bEnd)
{
    return aStart <= bEnd && bStart <= aEnd;
}

}

WindowStatuses buildWindowStatuses(const GcProfileData &gcData,
                                   const std::vector<ExcludedRegion> &exclusions,
                                   const DiploidRegions &diploidRegions)
{
    const bool checkDiploid = !diploidRegions.empty();

    // ---- SuppliedExcludedRegions 建構子：依染色體分組後排序 ----
    std::vector<std::string> exKeys;
    std::vector<std::vector<ExcludedRegion>> exLists;
    for(const ExcludedRegion &r : exclusions)
    {
        const std::string shortName = lp::stripChrPrefix(r.chromosome);
        std::size_t idx = exKeys.size();
        for(std::size_t i = 0; i < exKeys.size(); ++i)
        {
            if(exKeys[i] == shortName){ idx = i; break; }
        }
        if(idx == exKeys.size()){ exKeys.push_back(shortName); exLists.emplace_back(); }
        exLists[idx].push_back(r);
    }
    for(auto &v : exLists)
    {
        std::stable_sort(v.begin(), v.end(), [](const ExcludedRegion &a, const ExcludedRegion &b){
            if(a.start != b.start){ return a.start < b.start; }
            return a.end < b.end;
        });
    }

    // ---- findIntersections：帶狀態的單向掃描（見 behaviour-contract.md）----
    // 以 (染色體短名, gcProfile.start) 代表被排除的 window；每條染色體的 start 唯一，
    // 與 Java 端 GCProfile 的值相等判定等價。
    std::set<std::pair<std::string, int>> toExclude;
    for(std::size_t ci = 0; ci < gcData.chromosomes.size(); ++ci)
    {
        const std::string &shortName = gcData.chromosomes[ci];
        std::size_t ei = exKeys.size();
        for(std::size_t i = 0; i < exKeys.size(); ++i)
        {
            if(exKeys[i] == shortName){ ei = i; break; }
        }
        if(ei == exKeys.size()){ continue; }               // filterRegions == null

        const std::vector<GcProfile> &profiles = gcData.byChromosome[ci];
        long indexOfLastAffectedProfile = -1;
        for(const ExcludedRegion &region : exLists[ei])
        {
            bool keepSearching = true;
            bool found = false;
            for(std::size_t j = static_cast<std::size_t>(indexOfLastAffectedProfile + 1);
                j < profiles.size() && keepSearching; ++j)
            {
                const GcProfile &p = profiles[j];
                if(overlaps(region.start, region.end, p.start, p.end))
                {
                    toExclude.insert({shortName, p.start});
                    indexOfLastAffectedProfile = static_cast<long>(j);
                    found = true;
                }
                else
                {
                    keepSearching = !found;
                }
            }
        }
    }

    // ---- 逐 window 建 WindowStatus ----
    WindowStatuses out;
    out.chromosomes = gcData.chromosomes;
    out.byChromosome.resize(gcData.chromosomes.size());
    for(std::size_t ci = 0; ci < gcData.chromosomes.size(); ++ci)
    {
        const std::string &shortName = gcData.chromosomes[ci];
        const std::vector<DiploidStatus> *ds = diploidRegions.find(shortName);
        const std::vector<GcProfile> &profiles = gcData.byChromosome[ci];
        std::vector<WindowStatus> &dst = out.byChromosome[ci];
        dst.reserve(profiles.size());
        for(const GcProfile &g : profiles)
        {
            WindowStatus w;
            w.chromosome = g.chromosome;
            w.start = g.start;
            w.end = g.end;
            w.excluded = toExclude.count({shortName, g.start}) > 0;
            w.unmappable = !g.isMappable();
            w.nonDiploid = false;
            if(checkDiploid)
            {
                const std::size_t index = static_cast<std::size_t>(indexFor(g.start));
                const std::size_t size = (ds == nullptr) ? 0 : ds->size();
                w.nonDiploid = (index < size) ? !(*ds)[index].isDiploid : true;
            }
            dst.push_back(w);
        }
    }
    return out;
}

}
