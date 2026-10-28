#include "PurpleScaleCheck.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#include "../common/HumanChromosome.h"

namespace purple {
namespace {

constexpr double PI = 3.14159265358979323846;

// 與 PurpleWriters.cpp 的 number() 相同的 "%.4f"（含 0.00005 的特例），再讀回 double。
double round4(double value){
    if(value == 0.00005){ return 0.0; }
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);
    return std::strtod(buffer, nullptr);
}

// 與 PurpleWriters.cpp 的 minor()／major() 相同（cnv.somatic.tsv 的兩欄）。
double minorAllele(const PurpleCopyNumber &c){
    if(c.averageActualBaf < 0.5 - 1e-10){ return 0; }
    return std::max(0.0, (1 - c.averageActualBaf) * c.averageTumorCopyNumber);
}

bool autosome(const std::string &chromosome){
    const std::string name = lp::stripChrPrefix(chromosome);
    if(name.empty() || name.find_first_not_of("0123456789") != std::string::npos){ return false; }
    const int n = std::atoi(name.c_str());
    return n >= 1 && n <= 22;
}

double tdf(double purity, double ploidy){ return purity * ploidy / (purity * ploidy + 2 * (1 - purity)); }

std::string fixed(double value, int digits){
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.*f", digits, value);
    return buffer;
}

}

ScaleCheck checkScale(const std::vector<PurpleCopyNumber> &copyNumbers, const FittedPurity &selected,
        const std::vector<FittedPurity> &fits, bool f0Valid, double f0){
    ScaleCheck check;
    check.ploidy = round4(selected.ploidy);
    const double psi = check.ploidy;
    std::vector<double> yMajor, yMinor, weight;
    for(const auto &c : copyNumbers){
        if(!autosome(c.chromosome) || c.end - c.start < 1000000 || c.bafCount < 20){ continue; }
        const double minor = round4(minorAllele(c));
        const double major = round4(c.averageTumorCopyNumber - minorAllele(c));
        yMajor.push_back(major / psi);
        yMinor.push_back(minor / psi);
        weight.push_back(static_cast<double>(c.end - c.start));
    }
    check.segments = static_cast<int>(weight.size());
    if(weight.empty() || psi <= 0){ return check; }
    double totalWeight = 0;
    for(double w : weight){ totalWeight += w; }
    const auto amplitude = [&](const std::vector<double> &y, double s){
        double re = 0, im = 0;
        for(size_t i = 0; i < y.size(); ++i){
            const double angle = 2 * PI * s * y[i];
            re += weight[i] * std::cos(angle);
            im += weight[i] * std::sin(angle);
        }
        return std::sqrt(re * re + im * im) / totalWeight;
    };
    for(int k = 150; k <= 1000; ++k){
        const double s = k / 100.0;
        const double a = 0.5 * (amplitude(yMajor, s) + amplitude(yMinor, s));
        if(k == 150 || a > check.aMax){ check.aMax = a; check.sMax = s; }
    }
    check.computed = true;
    if(std::abs(check.sMax / (2 * psi) - 1) <= 0.25){ check.decision = "double"; }
    else if(std::abs(check.sMax / psi - 1) <= 0.25){ check.decision = "keep"; }
    else{ check.decision = "unclear"; }
    if(check.decision == "double" && psi >= 2 && f0Valid){
        for(const auto &fit : fits){
            if(fit.ploidy < 1.75 * psi || fit.ploidy > 2.25 * psi){ continue; }
            if(std::abs(tdf(fit.purity, fit.ploidy) - f0) > 0.05 + 1e-12){ continue; }
            if(!check.hasCandidate || fit.score < check.candidate.score){ check.candidate = fit; check.hasCandidate = true; }
        }
        check.applied = check.hasCandidate;
    }
    return check;
}

void writeScaleCheckTsv(const std::string &outputDir, const std::string &sampleId, const ScaleCheck &check,
        const FittedPurity &selected){
    if(outputDir.empty()){ return; }
    const std::string path = outputDir + (outputDir.back() == '/' ? "" : "/") + sampleId + ".purple.scale_check.tsv";
    std::ofstream out(path);
    if(!out){ throw std::runtime_error("unable to open output: " + path); }
    out << "segments\tploidy\tsMax\taMax\tdecision\tapplied\tselectedPurity\tselectedPloidy\tcandidatePurity\tcandidatePloidy\tcandidateScore\n";
    out << check.segments << '\t' << fixed(check.ploidy, 4) << '\t' << fixed(check.sMax, 2) << '\t' << fixed(check.aMax, 4) << '\t'
        << (check.computed ? check.decision : "NA") << '\t' << (check.applied ? "true" : "false") << '\t'
        << fixed(selected.purity, 4) << '\t' << fixed(selected.ploidy, 4) << '\t'
        << (check.hasCandidate ? fixed(check.candidate.purity, 4) : "NA") << '\t'
        << (check.hasCandidate ? fixed(check.candidate.ploidy, 4) : "NA") << '\t'
        << (check.hasCandidate ? fixed(check.candidate.score, 4) : "NA") << '\n';
}

}
