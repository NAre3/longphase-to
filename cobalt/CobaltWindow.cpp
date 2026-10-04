#include "CobaltWindow.h"

#include "CobaltConstants.h"
#include "../common/HumanChromosome.h"

namespace cobalt {

namespace {

// HumanChromosome.isAutosome()：1..22 為真，X／Y 為假。
// 實作於 cobalt/ 而非 common/，因為這是 COBALT 專屬需求；
// 動 common/ 會使 AMBER 的 check_no_regression.sh 需要重跑。
bool isAutosome(const std::string &shortName)
{
    if(shortName == "X" || shortName == "Y"){ return false; }
    for(char c : shortName){ if(c < '0' || c > '9'){ return false; } }
    return !shortName.empty();
}

inline int indexFor(int position){ return position / WINDOW_SIZE; }

}

std::vector<CobaltWindow> buildWindows(const std::vector<DepthReading> &depths,
                                       const WindowStatuses &statuses,
                                       const GcProfileData &gcData,
                                       GcPailsList &pails,
                                       WindowBuildStats &stats)
{
    std::vector<CobaltWindow> out;
    out.reserve(depths.size());

    for(const DepthReading &d : depths)
    {
        const std::string shortName = lp::stripChrPrefix(d.chromosome);

        // --- WindowStatuses.exclude ---
        bool isExcluded = false;
        {
            const std::vector<WindowStatus> *st = nullptr;
            for(std::size_t i = 0; i < statuses.chromosomes.size(); ++i)
            {
                if(statuses.chromosomes[i] == shortName){ st = &statuses.byChromosome[i]; break; }
            }
            const std::size_t idx = static_cast<std::size_t>(indexFor(d.startPosition));
            if(st != nullptr && idx < st->size()){ isExcluded = (*st)[idx].maskedOut(); }
        }
        const bool isInTargetRegion = true;                 // WholeGenome：Scope.onTarget 恆真（G15）

        CobaltWindow w;
        w.chromosomeShort = shortName;
        w.position = d.startPosition;
        w.readDepth = d.readDepth;
        w.gcContent = d.readGcContent;
        w.isInExcludedRegion = isExcluded;
        w.isInTargetRegion = isInTargetRegion;
        w.gcBucket = -1;
        w.gcReplaced = false;

        // --- correctedByReferenceValue ---
        const bool include = !w.isInExcludedRegion && w.isInTargetRegion;
        if(!(w.readDepth >= 1.0 || w.readDepth < 0) && include)
        {
            const std::vector<GcProfile> *profiles = gcData.find(shortName);
            const std::size_t idx = static_cast<std::size_t>(indexFor(w.position));
            if(profiles != nullptr && idx < profiles->size())
            {
                w.gcContent = (*profiles)[idx].gcContent;
                w.gcReplaced = true;
                ++stats.gcReplacedCount;
            }
            else
            {
                ++stats.referenceLookupFailures;            // 對應 catch(Exception) → 保留原值
            }
        }

        // --- bucketed ---
        GcPail &bucket = pails.getGcPail(w.gcContent);       // 對所有 window 都呼叫
        if(!w.isInExcludedRegion)
        {
            if(isAutosome(shortName) && w.isInTargetRegion && w.readDepth > 0)
            {
                bucket.addReading(w.readDepth);
                ++stats.bucketedReadings;
            }
            w.gcBucket = bucket.gc();
        }
        // excluded 者 gcBucket 維持 -1（Java 端的 null）

        out.push_back(std::move(w));
    }
    return out;
}

}
