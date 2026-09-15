#include "PurpleSegmentation.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

struct Cluster {
    std::string chromosome;
    int start = 0;
    int end = 0;
    std::vector<PcfPosition> positions;
};

int chromosomeRank(const std::string &chromosome){
    const std::string value = lp::stripChrPrefix(chromosome);
    if(value == "X"){ return 23; }
    if(value == "Y"){ return 24; }
    try { return std::stoi(value); } catch(...) { return 1000; }
}

std::unordered_map<std::string, int> readLengths(const std::string &referenceFasta){
    std::ifstream input(referenceFasta + ".fai");
    if(!input){ throw std::runtime_error("unable to open reference index: " + referenceFasta + ".fai"); }
    std::unordered_map<std::string, int> result;
    std::string chromosome;
    int length = 0;
    while(input >> chromosome >> length){
        std::string rest;
        std::getline(input, rest);
        if(lp::isHumanChromosome(chromosome)){ result[chromosome] = length; }
    }
    return result;
}

int centromere(const std::string &chromosome){
    static const std::unordered_map<std::string, int> values = {
        {"1",123605523},{"2",93139351},{"3",92214016},{"4",50726026},{"5",48272854},{"6",59191911},
        {"7",59498944},{"8",44955505},{"9",44377363},{"10",40640102},{"11",52751711},{"12",35977330},
        {"13",17025624},{"14",17086762},{"15",18362627},{"16",37295920},{"17",24849830},{"18",18161053},
        {"19",25844927},{"20",28237290},{"21",11890184},{"22",14004553},{"X",60509061},{"Y",10430492}
    };
    return values.at(lp::stripChrPrefix(chromosome));
}

std::vector<PcfPosition> mergePcf(const std::vector<PcfPosition> &ratio, const std::vector<PcfPosition> &baf){
    std::vector<PcfPosition> result;
    result.reserve(ratio.size() + baf.size());
    std::size_t i = 0, j = 0;
    while(i < ratio.size() || j < baf.size()){
        if(i == ratio.size()){ result.push_back(baf[j++]); }
        else if(j == baf.size() || ratio[i].position < baf[j].position){ result.push_back(ratio[i++]); }
        else if(baf[j].position < ratio[i].position){ result.push_back(baf[j++]); }
        else { result.push_back(ratio[i++]); ++j; }
    }
    return result;
}

int firstValidRatio(int position, int index, const std::vector<const CobaltRatio *> &ratios){
    const int windowStart = ((position - 1) / 1000) * 1000 + 1;
    const int minimum = windowStart - 999;
    for(int i = index; i >= 0; --i){
        if(ratios[static_cast<std::size_t>(i)]->position <= minimum && ratios[static_cast<std::size_t>(i)]->tumorGcRatio > -1){
            return ratios[static_cast<std::size_t>(i)]->position + 1;
        }
    }
    return minimum;
}

std::vector<Cluster> clustersForChromosome(const std::string &chromosome, const std::vector<PcfPosition> &positions,
                                            const std::vector<const CobaltRatio *> &ratios){
    std::vector<Cluster> result;
    int cobaltIndex = 0;
    Cluster *cluster = nullptr;
    bool lastWasAmberEnd = false;
    for(std::size_t i = 0; i < positions.size(); ++i){
        const PcfPosition &position = positions[i];
        if(position.position == 1){ continue; }
        while(cobaltIndex < static_cast<int>(ratios.size()) - 1 && ratios[static_cast<std::size_t>(cobaltIndex)]->position < position.position){
            ++cobaltIndex;
        }
        const int earliest = firstValidRatio(position.position, cobaltIndex, ratios);
        bool canSegment = true;
        if(cluster != nullptr){
            if(earliest > cluster->end){
                if(position.source == PcfSource::TUMOR_BAF){
                    const bool useAmber = lastWasAmberEnd && position.isSegmentStart() && i < positions.size() - 2;
                    if(!useAmber){ continue; }
                }
            }else{
                canSegment = false;
            }
        }else{
            canSegment = position.source != PcfSource::TUMOR_BAF;
        }
        if(cluster == nullptr || canSegment){
            result.push_back(Cluster{chromosome, earliest, 0, {}});
            cluster = &result.back();
        }
        cluster->end = position.position;
        lastWasAmberEnd = false;
        cluster->positions.push_back(position);
        if(position.source == PcfSource::TUMOR_BAF && position.isSegmentEnd()){ lastWasAmberEnd = true; }
    }
    return result;
}

std::vector<SupportSegment> segmentsForChromosome(const std::string &chromosome, int length, const std::vector<Cluster> &clusters){
    std::vector<SupportSegment> result;
    SupportSegment segment{chromosome, 1, 0, true, SegmentSupport::TELOMERE, false, 1, 1};
    for(const Cluster &cluster : clusters){
        if(cluster.positions.empty()){ continue; }
        const PcfPosition &first = cluster.positions.front();
        segment.end = first.position - 1;
        result.push_back(segment);
        int minStart = first.position;
        int maxStart = first.position;
        for(const PcfPosition &position : cluster.positions){
            if(position.source == PcfSource::TUMOR_RATIO || position.source == PcfSource::REFERENCE_RATIO){
                minStart = std::min(minStart, position.minPosition);
                maxStart = std::max(maxStart, position.maxPosition);
            }
        }
        // PurpleSupportSegmentFactory.createPcfSegment always marks PCF-only
        // segments as ratio-supported, even when the cluster contains only BAF PCFs.
        segment = SupportSegment{chromosome, first.position, 0, true, SegmentSupport::NONE, false, minStart, maxStart};
    }
    segment.end = length;
    result.push_back(segment);
    for(std::size_t i = 1; i < result.size(); ++i){ result[i].minStart = std::max(result[i].minStart, result[i - 1].end + 1); }

    const int center = centromere(chromosome);
    std::vector<SupportSegment> withCentromere;
    withCentromere.reserve(result.size() + 1);
    for(const SupportSegment &original : result){
        if(center >= original.start && center <= original.end){
            if(original.start == center){
                SupportSegment centered = original;
                centered.support = SegmentSupport::CENTROMERE;
                withCentromere.push_back(centered);
            }else{
                SupportSegment left = original;
                left.end = center - 1;
                left.maxStart = std::min(left.maxStart, left.end);
                SupportSegment right = original;
                right.start = center;
                right.minStart = center;
                right.maxStart = center;
                right.support = SegmentSupport::CENTROMERE;
                withCentromere.push_back(left);
                withCentromere.push_back(right);
            }
        }else{
            withCentromere.push_back(original);
        }
    }
    return withCentromere;
}

}

std::vector<SupportSegment> createSupportSegments(const InputData &data, const std::string &referenceFasta){
    const auto lengths = readLengths(referenceFasta);
    std::unordered_map<std::string, std::vector<PcfPosition>> amberByChr, cobaltByChr;
    std::unordered_map<std::string, std::vector<const CobaltRatio *>> ratiosByChr;
    for(const auto &position : data.amberPcf){ amberByChr[position.chromosome].push_back(position); }
    for(const auto &position : data.cobaltTumorPcf){ cobaltByChr[position.chromosome].push_back(position); }
    for(const auto &ratio : data.ratios){ ratiosByChr[ratio.chromosome].push_back(&ratio); }

    std::vector<std::string> chromosomes;
    for(const auto &entry : cobaltByChr){ chromosomes.push_back(entry.first); }
    for(const auto &entry : amberByChr){ if(cobaltByChr.count(entry.first) == 0){ chromosomes.push_back(entry.first); } }
    std::sort(chromosomes.begin(), chromosomes.end(), [](const std::string &a, const std::string &b){ return chromosomeRank(a) < chromosomeRank(b); });

    std::vector<SupportSegment> result;
    for(const std::string &chromosome : chromosomes){
        const auto positions = mergePcf(cobaltByChr[chromosome], amberByChr[chromosome]);
        const auto clusters = clustersForChromosome(chromosome, positions, ratiosByChr[chromosome]);
        const auto chromosomeSegments = segmentsForChromosome(chromosome, lengths.at(chromosome), clusters);
        result.insert(result.end(), chromosomeSegments.begin(), chromosomeSegments.end());
    }
    return result;
}

void dumpSupportSegments(const std::vector<SupportSegment> &segments){
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(segments.size());
    for(const auto &x : segments){
        rows.push_back({x.chromosome, x.start, x.chromosome + "\t" + std::to_string(x.start) + "\t" + std::to_string(x.end) + "\t" +
                (x.ratioSupport ? "true" : "false") + "\t" + segmentSupportName(x.support) + "\t" +
                (x.svCluster ? "true" : "false") + "\t" + std::to_string(x.minStart) + "\t" + std::to_string(x.maxStart)});
    }
    lp::CpDump::write("CP-P2-support-segments", "chromosome\tstart\tend\tratioSupport\tsupport\tsvCluster\tminStart\tmaxStart", rows);
}

}
