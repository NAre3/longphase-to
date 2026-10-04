#include "PurpleLohFill.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

constexpr double EPS = 1e-10;
constexpr double MIN_COVERAGE = 0.9;
constexpr double MAX_NON_LOH_COVERAGE = 0.1;
constexpr double MAX_RETENTION = 0.3;

bool greaterOrEqual(double a, double b){ return a - b > -EPS; }

bool isAutosome(const std::string &stripped){
    if(stripped.empty() || stripped.size() > 2){ return false; }
    for(char c : stripped){ if(c < '0' || c > '9'){ return false; } }
    const int n = std::stoi(stripped);
    return n >= 1 && n <= 22;
}

// 同 PurpleFitting.cpp 的 selectFittingRegions，但不要求 bafCount > 0，且限體染色體。
bool fillable(const ObservedRegion &region){
    return isAutosome(lp::stripChrPrefix(region.segment.chromosome)) &&
           region.germlineStatus == GermlineStatus::DIPLOID &&
           greaterOrEqual(region.observedTumorRatio, 0) && !(region.observedTumorRatio - 3 > EPS);
}

double coverage(const ObservedRegion &region, const LohFillInput &input){
    const auto it = input.lohSpans.find(lp::stripChrPrefix(region.segment.chromosome));
    if(it == input.lohSpans.end()){ return 0.0; }
    const int start = region.segment.start, end = region.segment.end;
    long long overlap = 0;
    for(const auto &span : it->second){
        const int lo = std::max(span.first, start), hi = std::min(span.second, end);
        if(hi >= lo){ overlap += static_cast<long long>(hi) - lo + 1; }
    }
    return static_cast<double>(overlap) / (static_cast<double>(end) - start + 1);
}

int panelCount(const ObservedRegion &region, const LohFillInput &input){
    const auto it = input.panelPositions.find(lp::stripChrPrefix(region.segment.chromosome));
    if(it == input.panelPositions.end()){ return 0; }
    const auto &v = it->second;
    return static_cast<int>(std::upper_bound(v.begin(), v.end(), region.segment.end) -
                            std::lower_bound(v.begin(), v.end(), region.segment.start));
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

std::map<std::string, std::vector<int>> loadPanelPositions(const std::string &path){
    gzFile file = gzopen(path.c_str(), "rb");
    if(file == nullptr){ throw std::runtime_error("unable to open germline site panel: " + path); }
    std::map<std::string, std::vector<int>> result;
    char buffer[4096];
    bool header = true;
    while(gzgets(file, buffer, sizeof(buffer)) != nullptr){
        if(header){ header = false; continue; }
        const char *tab = std::strchr(buffer, '\t');
        if(tab == nullptr){ continue; }
        const std::string chromosome = lp::stripChrPrefix(std::string(buffer, tab - buffer));
        if(!isAutosome(chromosome)){ continue; }
        result[chromosome].push_back(std::atoi(tab + 1));
    }
    gzclose(file);
    for(auto &entry : result){ std::sort(entry.second.begin(), entry.second.end()); }
    return result;
}

std::map<std::string, std::vector<std::pair<int, int>>> loadLohBed(const std::string &path){
    std::ifstream in(path);
    if(!in){ throw std::runtime_error("unable to open LOH bed: " + path); }
    std::map<std::string, std::vector<std::pair<int, int>>> result;
    std::string line;
    while(std::getline(in, line)){
        std::istringstream fields(line);
        std::string chromosome;
        int start = 0, end = 0;
        if(!(fields >> chromosome >> start >> end)){ continue; }
        chromosome = lp::stripChrPrefix(chromosome);
        if(isAutosome(chromosome)){ result[chromosome].push_back({start, end}); }
    }
    return result;
}

LohFillResult applyLohFill(std::vector<ObservedRegion> &regions, const LohFillInput &input){
    LohFillResult result;
    std::vector<double> cov(regions.size(), 0.0);
    std::vector<int> pan(regions.size(), 0);
    long long sumBaf = 0, sumPanel = 0;
    for(std::size_t i = 0; i < regions.size(); ++i){
        const auto &region = regions[i];
        if(!isAutosome(lp::stripChrPrefix(region.segment.chromosome))){ continue; }
        cov[i] = coverage(region, input);
        pan[i] = panelCount(region, input);
        if(fillable(region) && cov[i] < MAX_NON_LOH_COVERAGE && region.bafCount > 0){
            sumBaf += region.bafCount;
            sumPanel += pan[i];
            ++result.retSegments;
        }
    }
    if(sumPanel <= 0){
        std::cerr << "PURPLE LOH fill: no non-LOH reference segment, skipped\n";
        return result;
    }
    result.ret = static_cast<double>(sumBaf) / static_cast<double>(sumPanel);
    for(std::size_t i = 0; i < regions.size(); ++i){
        auto &region = regions[i];
        if(!fillable(region) || !(cov[i] >= MIN_COVERAGE) || pan[i] <= 0){ continue; }
        const double expected = result.ret * pan[i];
        const double retention = region.bafCount / expected;
        if(!(retention < MAX_RETENTION)){ continue; }
        LohFillRecord record;
        record.chromosome = region.segment.chromosome;
        record.start = region.segment.start;
        record.end = region.segment.end;
        record.originalBafCount = region.bafCount;
        record.originalObservedBaf = region.observedBaf;
        record.panelCount = pan[i];
        record.coverage = cov[i];
        record.expected = expected;
        record.retention = retention;
        // Python round() 為四捨六入五成雙；nearbyint 在預設捨入模式下相同。
        record.filledBafCount = std::max(1, static_cast<int>(std::nearbyint(expected)));
        region.bafCount = record.filledBafCount;
        region.observedBaf = 1.0;
        result.filled.push_back(record);
    }
    return result;
}

void writeLohFillTsv(const std::string &outputDir, const std::string &sampleId, const LohFillResult &result){
    if(outputDir.empty()){ return; }
    const std::string path = joinPath(outputDir, sampleId + ".purple.loh_fill.tsv");
    std::ofstream out(path);
    if(!out){ throw std::runtime_error("unable to open output: " + path); }
    out << "#retention=" << lp::CpDump::num(result.ret) << "\tretentionSegments=" << result.retSegments
        << "\tfilledSegments=" << result.filled.size() << '\n';
    out << "chromosome\tstart\tend\toriginalBafCount\toriginalObservedBAF\tpanelSites\tlohCoverage\texpectedBafCount\tretention\tfilledBafCount\n";
    for(const auto &r : result.filled){
        out << r.chromosome << '\t' << r.start << '\t' << r.end << '\t' << r.originalBafCount << '\t'
            << fixed(r.originalObservedBaf, 4) << '\t' << r.panelCount << '\t' << fixed(r.coverage, 4) << '\t'
            << fixed(r.expected, 2) << '\t' << fixed(r.retention, 4) << '\t' << r.filledBafCount << '\n';
    }
}

void dumpLohFill(const LohFillResult &result){
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(result.filled.size());
    for(const auto &r : result.filled){
        rows.push_back({r.chromosome, r.start, r.chromosome + "\t" + std::to_string(r.start) + "\t" + std::to_string(r.end) + "\t" +
            std::to_string(r.originalBafCount) + "\t" + lp::CpDump::num(r.originalObservedBaf) + "\t" + std::to_string(r.panelCount) + "\t" +
            lp::CpDump::num(r.coverage) + "\t" + lp::CpDump::num(r.expected) + "\t" + lp::CpDump::num(r.retention) + "\t" +
            std::to_string(r.filledBafCount) + "\t" + lp::CpDump::num(result.ret)});
    }
    lp::CpDump::write("CP-P3b-loh-fill", "chromosome\tstart\tend\toriginalBafCount\toriginalObservedBAF\tpanelSites\tlohCoverage\texpectedBafCount\tretention\tfilledBafCount\tret", rows);
}

}
