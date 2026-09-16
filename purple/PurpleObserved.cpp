#include "PurpleObserved.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

struct Interval { const char *chromosome; int start; int end; };

const Interval IMMUNE_REGIONS[] = {
    {"chr2",88857161,90315836},{"chr6",29937532,29950870},{"chr6",31263749,31277092},
    {"chr6",31348872,31372067},{"chr7",38239580,38367882},{"chr7",142299177,142812869},
    {"chr14",21622293,22552156},{"chr14",105586437,106879844},{"chr22",22026076,22922913}
};

int centromere(const std::string &chromosome){
    static const std::unordered_map<std::string, int> values = {
        {"1",123605523},{"2",93139351},{"3",92214016},{"4",50726026},{"5",48272854},{"6",59191911},
        {"7",59498944},{"8",44955505},{"9",44377363},{"10",40640102},{"11",52751711},{"12",35977330},
        {"13",17025624},{"14",17086762},{"15",18362627},{"16",37295920},{"17",24849830},{"18",18161053},
        {"19",25844927},{"20",28237290},{"21",11890184},{"22",14004553},{"X",60509061},{"Y",10430492}
    };
    return values.at(lp::stripChrPrefix(chromosome));
}

bool within(int innerStart, int innerEnd, int outerStart, int outerEnd){
    return innerStart <= innerEnd && innerStart >= outerStart && innerEnd <= outerEnd;
}

std::vector<SupportSegment> split(const SupportSegment &segment, int position){
    if(position <= segment.start || position > segment.end){ return {segment}; }
    SupportSegment left = segment;
    left.end = position - 1;
    left.maxStart = std::min(left.maxStart, left.end);
    SupportSegment right = segment;
    right.start = position;
    right.support = SegmentSupport::EXCL;
    right.minStart = position;
    right.maxStart = position;
    return {left, right};
}

std::vector<SupportSegment> refineExcluded(const std::vector<SupportSegment> &segments){
    std::vector<SupportSegment> current = segments;
    for(const Interval &excluded : IMMUNE_REGIONS){
        std::vector<SupportSegment> next;
        for(const SupportSegment &segment : current){
            if(segment.chromosome != excluded.chromosome){ next.push_back(segment); continue; }
            for(const SupportSegment &first : split(segment, excluded.start)){
                const auto second = split(first, excluded.end + 1);
                next.insert(next.end(), second.begin(), second.end());
            }
        }
        current.swap(next);
    }
    std::sort(current.begin(), current.end(), [](const SupportSegment &a, const SupportSegment &b){
        if(a.chromosome != b.chromosome){
            const std::string ac = lp::stripChrPrefix(a.chromosome), bc = lp::stripChrPrefix(b.chromosome);
            const int ar = ac == "X" ? 23 : ac == "Y" ? 24 : std::stoi(ac);
            const int br = bc == "X" ? 23 : bc == "Y" ? 24 : std::stoi(bc);
            return ar < br;
        }
        return a.start < b.start;
    });
    return current;
}

double median(std::vector<double> values){
    if(values.empty()){ return 0; }
    std::sort(values.begin(), values.end());
    const std::size_t middle = values.size() / 2;
    return values.size() % 2 == 0 ? (values[middle - 1] + values[middle]) * 0.5 : values[middle];
}

bool excludedStatus(const SupportSegment &segment){
    for(const Interval &region : IMMUNE_REGIONS){
        if(segment.chromosome == region.chromosome && within(segment.start, segment.end, region.start, region.end)){ return true; }
    }
    return false;
}

bool diploidBafChromosome(const std::string &chromosome, Gender gender){
    const std::string value = lp::stripChrPrefix(chromosome);
    if(value == "Y"){ return false; }
    if(value == "X"){ return gender == Gender::FEMALE; }
    return lp::isHumanChromosome(chromosome);
}

GermlineStatus germlineStatus(const SupportSegment &segment, Gender gender,
        double normalRatio, double tumorRatio, int depthWindowCount){
    if(excludedStatus(segment)){ return GermlineStatus::EXCLUDED; }
    if(within(segment.start, segment.end, centromere(segment.chromosome) - 2000000,
              centromere(segment.chromosome) + 2000000)){ return GermlineStatus::CENTROMETIC; }
    const std::string chromosome = lp::stripChrPrefix(segment.chromosome);
    if(depthWindowCount == 0 || (chromosome == "Y" && gender == Gender::FEMALE)){ return GermlineStatus::UNKNOWN; }
    const double adjustment = (chromosome == "X" || chromosome == "Y") && gender == Gender::MALE ? 0.5 : 1.0;
    if(normalRatio < 0.1 * adjustment - 1e-10 && tumorRatio < 0.1 * adjustment - 1e-10){ return GermlineStatus::HOM_DELETION; }
    if(normalRatio < 0.7 * adjustment - 1e-10){ return GermlineStatus::HET_DELETION; }
    if(normalRatio < 0.85 * adjustment - 1e-10){ return GermlineStatus::LIKELY_DIPLOID; }
    if(normalRatio < 1.15 * adjustment - 1e-10){ return GermlineStatus::DIPLOID; }
    if(normalRatio < 1.3 * adjustment - 1e-10){ return GermlineStatus::LIKELY_DIPLOID; }
    if(normalRatio < 2.2 * adjustment - 1e-10){ return GermlineStatus::AMPLIFICATION; }
    return GermlineStatus::NOISE;
}

}

std::vector<ObservedRegion> createObservedRegions(const InputData &data, const std::vector<SupportSegment> &segments){
    const auto refined = refineExcluded(segments);
    std::unordered_map<std::string, std::vector<const AmberBaf *>> bafByChr;
    std::unordered_map<std::string, std::vector<const CobaltRatio *>> ratioByChr;
    for(const auto &baf : data.bafs){ bafByChr[baf.chromosome].push_back(&baf); }
    for(const auto &ratio : data.ratios){ ratioByChr[ratio.chromosome].push_back(&ratio); }
    std::unordered_map<std::string, std::size_t> bafIndex, ratioIndex;

    std::vector<ObservedRegion> result;
    result.reserve(refined.size());
    for(const SupportSegment &segment : refined){
        std::vector<double> bafValues;
        if(diploidBafChromosome(segment.chromosome, data.cobaltGender)){
            const auto &values = bafByChr[segment.chromosome];
            auto &index = bafIndex[segment.chromosome];
            while(index < values.size() && values[index]->position < segment.start){ ++index; }
            std::size_t cursor = index;
            while(cursor < values.size() && values[cursor]->position <= segment.end){
                // Match AmberBAF.tumorModifiedBAF() exactly.  The mathematically
                // equivalent max(x, 1-x) can differ by one ULP at .xxxx5 and
                // therefore change Java's four-decimal writer output.
                bafValues.push_back(0.5 + std::abs(values[cursor]->tumorBaf - 0.5));
                ++cursor;
            }
            index = cursor;
        }

        double referenceSum = 0, rawReferenceSum = 0, tumorContentSum = 0;
        int referenceCount = 0, tumorCount = 0;
        std::vector<double> tumorRatios;
        const auto &values = ratioByChr[segment.chromosome];
        auto &index = ratioIndex[segment.chromosome];
        while(index < values.size() && values[index]->position < segment.start){ ++index; }
        std::size_t cursor = index;
        while(cursor < values.size() && values[cursor]->position <= segment.end){
            const CobaltRatio &ratio = *values[cursor];
            const int windowEnd = ((ratio.position - 1) / 1000) * 1000 + 1000;
            if(windowEnd <= segment.end && ratio.referenceGcDiploidRatio >= 0){
                referenceSum += ratio.referenceGcDiploidRatio;
                rawReferenceSum += ratio.referenceGcRatio;
                ++referenceCount;
                if(ratio.tumorGcRatio > -1){
                    tumorRatios.push_back(ratio.tumorGcRatio);
                    tumorContentSum += ratio.tumorGcContent;
                    ++tumorCount;
                }
            }
            ++cursor;
        }
        index = cursor;

        ObservedRegion observed;
        observed.segment = segment;
        observed.bafCount = static_cast<int>(bafValues.size());
        observed.observedBaf = median(std::move(bafValues));
        observed.depthWindowCount = tumorCount;
        observed.observedTumorRatio = median(std::move(tumorRatios));
        observed.observedNormalRatio = referenceCount > 0 ? referenceSum / referenceCount : 0;
        observed.unnormalisedObservedNormalRatio = referenceCount > 0 ? rawReferenceSum / referenceCount : 0;
        observed.gcContent = tumorCount > 0 ? tumorContentSum / tumorCount : 0;

        observed.germlineStatus = germlineStatus(segment, data.cobaltGender, observed.observedNormalRatio,
                                                 observed.observedTumorRatio, tumorCount);
        result.push_back(std::move(observed));
    }

    for(std::size_t i = 0; i < result.size(); ++i){
        auto &target = result[i];
        if(target.segment.support != SegmentSupport::NONE || target.germlineStatus != GermlineStatus::DIPLOID){ continue; }
        for(std::size_t j = i; j-- > 0; ){
            const auto &prior = result[j];
            if(prior.germlineStatus == GermlineStatus::DIPLOID){ break; }
            target.segment.minStart = std::min(target.segment.minStart, prior.segment.start);
            if(prior.segment.support != SegmentSupport::NONE){ break; }
        }
    }
    return result;
}

void dumpObservedRegions(const std::vector<ObservedRegion> &regions){
    using lp::CpDump;
    std::vector<CpDump::Row> rows;
    rows.reserve(regions.size());
    for(const auto &x : regions){
        const auto &s = x.segment;
        rows.push_back({s.chromosome, s.start, s.chromosome + "\t" + std::to_string(s.start) + "\t" + std::to_string(s.end) + "\t" +
            (s.ratioSupport ? "true" : "false") + "\t" + segmentSupportName(s.support) + "\t" + std::to_string(x.bafCount) + "\t" +
            CpDump::num(x.observedBaf) + "\t" + std::to_string(x.depthWindowCount) + "\t" + CpDump::num(x.observedTumorRatio) + "\t" +
            CpDump::num(x.observedNormalRatio) + "\t" + CpDump::num(x.unnormalisedObservedNormalRatio) + "\t" + germlineStatusName(x.germlineStatus) +
            "\t" + (s.svCluster ? "true" : "false") + "\t" + CpDump::num(x.gcContent) + "\t" + std::to_string(s.minStart) + "\t" +
            std::to_string(s.maxStart) + "\t0\t0\t0\t0\t0\t0\t0\t0\t0"});
    }
    CpDump::write("CP-P3-observed-regions", "chromosome\tstart\tend\tratioSupport\tsupport\tbafCount\tobservedBAF\tdepthWindowCount\tobservedTumorRatio\tobservedNormalRatio\tunnormalisedObservedNormalRatio\tgermlineStatus\tsvCluster\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumberDeviation\tmajorAlleleCopyNumberDeviation\tdeviationPenalty\teventPenalty\trefNormalisedCopyNumber\ttumorCopyNumber\ttumorBAF\tfittedTumorCopyNumber\tfittedBAF", rows);
}

}
