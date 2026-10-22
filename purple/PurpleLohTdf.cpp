#include "PurpleLohTdf.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>

#include "../common/HumanChromosome.h"

namespace purple {
namespace {

constexpr double BIN_WIDTH = 0.005;
constexpr int SMOOTH_HALF_WIDTH = 3;
constexpr double PEAK_MIN_RELATIVE = 0.25;

bool isAutosome(const std::string &stripped){
    if(stripped.empty() || stripped.size() > 2){ return false; }
    for(char c : stripped){ if(c < '0' || c > '9'){ return false; } }
    const int n = std::stoi(stripped);
    return n >= 1 && n <= 22;
}

// 平滑值只在有資料的格上定義；鄰格無資料時視為 0。
double rightmostPeak(const std::vector<LohTdfSegment> &segments){
    std::map<int, long long> histogram;
    for(const auto &s : segments){ histogram[static_cast<int>(s.f / BIN_WIDTH)] += s.weight; }
    std::map<int, long long> smoothed;
    long long top = 0;
    for(const auto &entry : histogram){
        long long sum = 0;
        for(int d = -SMOOTH_HALF_WIDTH; d <= SMOOTH_HALF_WIDTH; ++d){
            const auto it = histogram.find(entry.first + d);
            if(it != histogram.end()){ sum += it->second; }
        }
        smoothed[entry.first] = sum;
        if(sum > top){ top = sum; }
    }
    auto valueAt = [&smoothed](int bin) -> long long {
        const auto it = smoothed.find(bin);
        return it == smoothed.end() ? 0 : it->second;
    };
    bool found = false;
    int best = 0;
    for(const auto &entry : smoothed){
        const int bin = entry.first;
        const long long value = entry.second;
        if(value >= PEAK_MIN_RELATIVE * top && value >= valueAt(bin - 1) && value >= valueAt(bin + 1)){
            if(!found || bin > best){ best = bin; found = true; }
        }
    }
    return (best + 0.5) * BIN_WIDTH;
}

std::string joinPath(const std::string &directory, const std::string &name){
    return directory + (directory.empty() || directory.back() == '/' ? "" : "/") + name;
}

std::string fixed(double value, int digits){
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}

}

double tumorDnaFraction(double purity, double ploidy){
    return purity * ploidy / (purity * ploidy + 2 * (1 - purity));
}

LohTdfEstimate estimateLohTdf(const std::vector<ObservedRegion *> &regions, const PhasedBafResult &phased){
    LohTdfEstimate estimate;
    double lengthSum = 0, weightedRatio = 0;
    for(const ObservedRegion *region : regions){
        if(!isAutosome(lp::stripChrPrefix(region->segment.chromosome))){ continue; }
        const double length = region->segment.end - region->segment.start + 1;
        lengthSum += length;
        weightedRatio += length * region->observedTumorRatio;
    }
    if(lengthSum <= 0){ return estimate; }
    estimate.meanRatio = weightedRatio / lengthSum;
    for(size_t i = 0; i < regions.size(); ++i){
        const ObservedRegion &region = *regions[i];
        if(!isAutosome(lp::stripChrPrefix(region.segment.chromosome))){ continue; }
        const bool filled = phased.lohFilled[i] != 0;
        if(!filled && phased.cellCount[i] <= 0){ continue; }
        LohTdfSegment s;
        s.chromosome = region.segment.chromosome;
        s.start = region.segment.start;
        s.end = region.segment.end;
        s.weight = region.bafCount;
        s.ratio = region.observedTumorRatio / estimate.meanRatio;
        s.baf = filled ? 1.0 : phased.phasedBaf[i];
        s.f = 1 - 2 * s.ratio * (1 - s.baf);
        s.filled = filled;
        if(s.f > 0 && s.f <= 1){ estimate.segments.push_back(s); }
    }
    if(estimate.segments.empty()){ return estimate; }
    estimate.f0 = rightmostPeak(estimate.segments);
    estimate.valid = true;
    return estimate;
}

std::vector<FittedPurity> restrictToTdfWindow(const std::vector<FittedPurity> &fits, double f0, double halfWidth){
    std::vector<FittedPurity> kept;
    for(const auto &fit : fits){
        if(std::abs(tumorDnaFraction(fit.purity, fit.ploidy) - f0) <= halfWidth){ kept.push_back(fit); }
    }
    return kept;
}

void writeLohTdfTsv(const std::string &outputDir, const std::string &sampleId, const LohTdfEstimate &estimate,
        double halfWidth, size_t candidatesKept, size_t candidatesTotal){
    if(outputDir.empty()){ return; }
    const std::string path = joinPath(outputDir, sampleId + ".purple.loh_tdf.tsv");
    std::ofstream out(path);
    if(!out){ throw std::runtime_error("unable to open output: " + path); }
    out << "#f0=" << (estimate.valid ? fixed(estimate.f0, 4) : "NA") << "\thalfWidth=" << fixed(halfWidth, 4)
        << "\tmeanRatio=" << fixed(estimate.meanRatio, 6) << "\tcandidatesKept=" << candidatesKept
        << "\tcandidatesTotal=" << candidatesTotal << '\n';
    out << "chromosome\tstart\tend\tbafCount\tratio\tbaf\tf\tlohFilled\n";
    for(const auto &s : estimate.segments){
        out << s.chromosome << '\t' << s.start << '\t' << s.end << '\t' << s.weight << '\t' << fixed(s.ratio, 6) << '\t'
            << fixed(s.baf, 6) << '\t' << fixed(s.f, 6) << '\t' << (s.filled ? 1 : 0) << '\n';
    }
}

}
