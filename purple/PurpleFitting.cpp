#include "PurpleFitting.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

constexpr double EPS = 1e-10;

double round2(double value){ return std::round(value * 100.0) / 100.0; }
bool lessThan(double a, double b){ return a - b < -EPS; }
bool lessOrEqual(double a, double b){ return a - b < EPS; }
bool greaterOrEqual(double a, double b){ return a - b > -EPS; }

double binomialCdf(int n, int k){
    double probability = std::ldexp(1.0, -n);
    double sum = probability;
    for(int i = 0; i < k; ++i){
        probability *= static_cast<double>(n - i) / static_cast<double>(i + 1);
        sum += probability;
    }
    return sum;
}

double expectedBaf(int averageDepth){
    const int minimum = averageDepth / 2;
    const int maximum = averageDepth * 3 / 2;
    double total = 0;
    for(int depth = minimum; depth < maximum; ++depth){
        int quantile = 0;
        while(quantile < depth && binomialCdf(depth, quantile) < 0.8){ ++quantile; }
        total += static_cast<double>(quantile) / depth;
    }
    return std::max(total / (maximum - minimum), 0.535);
}

double adjustedCopyNumber(double tumorRatio, double normalRatio, double purity, double normFactor){
    if(std::abs(tumorRatio) < EPS){ return 0; }
    return 2 * normalRatio + 2 * (tumorRatio - normalRatio * normFactor) / purity / normFactor;
}

double adjustedBafSimple(double copyNumber, double observed, double purity){
    if(lessOrEqual(copyNumber, 1)){ return 1; }
    const double totalObservations = purity * copyNumber + 2 * (1 - purity);
    return (observed * totalObservations - (1 - purity)) / purity / copyNumber;
}

double deviationProbability(double value, double ploidy){
    const double cdf = 0.5 * (1 + std::erf(value / std::sqrt(2.0)));
    return 2 * cdf - 1 + std::max(-0.5 - ploidy, 0.0);
}

double alleleDeviation(double purity, double normFactor, double ploidy){
    const double distance = lessThan(ploidy, -0.5) ? 0.5 : std::abs(ploidy - std::round(ploidy));
    const double deviationsPerPloidy = std::max(1.5, purity * normFactor / 2 / 0.05);
    return deviationProbability(distance * deviationsPerPloidy, ploidy);
}

double additionalPenalty(double minimum, double ploidy){
    return std::min(1.5, std::max(0.0, -1.5 * (ploidy - minimum)));
}

double majorDeviation(double purity, double normFactor, double ploidy){
    const double multiplier = ploidy > EPS && lessThan(ploidy, 1) ? std::max(1.0, 1.0 * (1 - ploidy)) : 1;
    return std::max(multiplier * alleleDeviation(purity, normFactor, ploidy) + additionalPenalty(1, ploidy), 0.1);
}

double minorDeviation(double purity, double normFactor, double ploidy){
    return std::max(alleleDeviation(purity, normFactor, ploidy) + additionalPenalty(0, ploidy), 0.1);
}

double estimateMinMaxBaf(double copyNumber, double minBaf, double maxBaf){
    const double majorMin = minBaf * copyNumber, majorMax = maxBaf * copyNumber;
    const double minorMin = copyNumber - majorMin, minorMax = copyNumber - majorMax;
    const double minorCeil = std::ceil(minorMax), majorCeil = std::ceil(majorMin);
    const auto sign = [](double x){ return x > 0 ? 1 : x < 0 ? -1 : 0; };
    const bool minorDiff = sign(minorCeil - minorMin) != sign(minorCeil - minorMax);
    const bool majorDiff = sign(majorCeil - majorMin) != sign(majorCeil - majorMax);
    if(!minorDiff && !majorDiff){ return -1; }
    if(!minorDiff){ return majorCeil / copyNumber; }
    if(!majorDiff){ return 1 - minorCeil / copyNumber; }
    const double half = copyNumber * 0.5;
    const double majorHigh = std::floor(majorMax) >= half ? std::floor(majorMax) : majorMax;
    double majorLow = majorHigh - 1;
    for(int minorInt = 0; minorInt < copyNumber; ++minorInt){
        const double major = copyNumber - minorInt;
        if(major > minorInt && major - minorInt <= 0.5){ majorLow = major; break; }
        if(major >= half){ majorLow = major; }
        else if(major < minorInt){
            const double previousFloor = std::floor(major + 1);
            if(previousFloor >= half){ majorLow = previousFloor; }
            break;
        }
    }
    majorLow = std::min(majorLow, majorHigh);
    const double diffHigh = majorHigh - (copyNumber - majorHigh);
    const double diffLow = majorLow - (copyNumber - majorLow);
    return std::max((diffHigh <= diffLow ? majorHigh : majorLow) / copyNumber, 0.5);
}

double impliedBaf(double copyNumber, double observed, double purity, double normFactor, double ambiguous){
    if(lessOrEqual(copyNumber, 0.1)){ return 1; }
    const auto adjusted = [&](double baf){ return adjustedBafSimple(copyNumber, baf, purity); };
    if(!lessOrEqual(observed, ambiguous)){ return adjusted(observed); }
    const double minBaf = std::clamp(adjusted(0.5), 0.0, 1.0);
    const double maxBaf = std::clamp(adjusted(observed), 0.0, 1.0);
    const double estimated = estimateMinMaxBaf(copyNumber, minBaf, maxBaf);
    if(estimated != -1){ return estimated; }
    const double majorMin = minBaf * copyNumber, majorMax = maxBaf * copyNumber;
    const double minDeviation = majorDeviation(purity, normFactor, majorMin) + minorDeviation(purity, normFactor, copyNumber - majorMin);
    const double maxDeviation = majorDeviation(purity, normFactor, majorMax) + minorDeviation(purity, normFactor, copyNumber - majorMax);
    return lessThan(minDeviation, maxDeviation) ? 0.5 : observed;
}

FittedPurity fitOne(double purity, double targetPloidy, double averageRatio, int totalBaf,
                    const std::vector<const ObservedRegion *> &regions, double ambiguous){
    const double normFactor = 2 * averageRatio / (2 - 2 * purity + targetPloidy * purity);
    double events = 0, deviations = 0, diploid = 0, averagePloidy = 0;
    for(const ObservedRegion *region : regions){
        const double copyNumber = adjustedCopyNumber(region->observedTumorRatio, 1, purity, normFactor);
        const double baf = impliedBaf(copyNumber, region->observedBaf, purity, normFactor, ambiguous);
        const double major = baf * copyNumber, minor = copyNumber - major;
        const double event = 1 + 0.4 * std::min(std::abs(major - 1) + std::abs(minor - 1),
                                               1 + std::abs(major - 2) + std::abs(minor - 2));
        const double deviation = (minorDeviation(purity, normFactor, minor) + majorDeviation(purity, normFactor, major)) * region->observedBaf;
        const double weight = static_cast<double>(region->bafCount) / totalBaf;
        events += event * weight;
        deviations += deviation * weight;
        averagePloidy += copyNumber * weight;
        if(greaterOrEqual(major, 0.8) && lessOrEqual(major, 1.2) && greaterOrEqual(minor, 0.8) && lessOrEqual(minor, 1.2)){
            diploid += weight;
        }
    }
    return FittedPurity{purity, normFactor, averagePloidy, events * deviations, diploid, 0};
}

std::vector<double> ploidies(){
    std::vector<double> values;
    for(double p = 1.0; lessThan(p, 3.0); p = round2(p + 0.02)){ values.push_back(p); }
    for(double p = 3.0; lessThan(p, 5.0); p = round2(p + 0.05)){ values.push_back(p); }
    for(double p = 5.0; lessThan(p, 8.0); p = round2(p + 0.1)){ values.push_back(p); }
    values.push_back(8.0);
    return values;
}

}

std::vector<const ObservedRegion *> selectFittingRegions(const std::vector<ObservedRegion> &regions){
    std::vector<const ObservedRegion *> result;
    for(const auto &region : regions){
        const std::string chromosome = lp::stripChrPrefix(region.segment.chromosome);
        if(region.bafCount <= 0 || !greaterOrEqual(region.observedTumorRatio, 0) || region.germlineStatus != GermlineStatus::DIPLOID ||
           region.observedTumorRatio - 3 > EPS || chromosome == "Y"){ continue; }
        result.push_back(&region);
    }
    return result;
}

std::vector<FittedPurity> fitPurityGrid(const std::vector<const ObservedRegion *> &regions, int averageTumorDepth, int threads){
    int totalBaf = 0;
    double weightedRatio = 0;
    for(const auto *region : regions){ totalBaf += region->bafCount; weightedRatio += region->bafCount * region->observedTumorRatio; }
    const double averageRatio = weightedRatio / totalBaf;
    const double ambiguous = expectedBaf(averageTumorDepth);
    std::vector<double> purities;
    for(double purity = 0.08; purity <= 1.0; purity = round2(purity + 0.01)){ purities.push_back(purity); }
    const auto ploidyValues = ploidies();
    std::vector<std::vector<FittedPurity>> byPurity(purities.size());
    std::atomic<std::size_t> next{0};
    const int workerCount = std::max(1, threads);
    std::vector<std::thread> workers;
    for(int worker = 0; worker < workerCount; ++worker){
        workers.emplace_back([&]{
            while(true){
                const std::size_t index = next.fetch_add(1);
                if(index >= purities.size()){ break; }
                auto &fits = byPurity[index];
                fits.reserve(ploidyValues.size());
                for(double ploidy : ploidyValues){ fits.push_back(fitOne(purities[index], ploidy, averageRatio, totalBaf, regions, ambiguous)); }
            }
        });
    }
    for(auto &worker : workers){ worker.join(); }
    std::vector<FittedPurity> result;
    result.reserve(purities.size() * ploidyValues.size());
    for(auto &fits : byPurity){ result.insert(result.end(), fits.begin(), fits.end()); }
    std::stable_sort(result.begin(), result.end(), [](const FittedPurity &a, const FittedPurity &b){ return a.score < b.score; });
    return result;
}

BestFit selectTumorOnlyBestFit(const std::vector<FittedPurity> &fits, const std::vector<ObservedRegion> &regions){
    const FittedPurity &lowest = fits.front();
    std::vector<const FittedPurity *> nearLowest;
    for(const auto &fit : fits){
        const double absolute = std::abs(fit.score - lowest.score);
        const double relative = std::abs(absolute / lowest.score);
        if(lessOrEqual(absolute, 0.0005) || lessOrEqual(relative, 0.1)){ nearLowest.push_back(&fit); }
    }
    FittedPurityScore score;
    score.minPurity = score.minPloidy = score.minDiploidProportion = std::numeric_limits<double>::max();
    for(const auto *fit : nearLowest){
        score.minPurity = std::min(score.minPurity, fit->purity); score.maxPurity = std::max(score.maxPurity, fit->purity);
        score.minPloidy = std::min(score.minPloidy, fit->ploidy); score.maxPloidy = std::max(score.maxPloidy, fit->ploidy);
        score.minDiploidProportion = std::min(score.minDiploidProportion, fit->diploidProportion);
        score.maxDiploidProportion = std::max(score.maxDiploidProportion, fit->diploidProportion);
    }

    int totalBaf = 0, highBaf = 0;
    for(const auto &region : regions){
        if(region.bafCount <= 0 || !greaterOrEqual(region.observedTumorRatio, 0) ||
           region.germlineStatus != GermlineStatus::DIPLOID || lp::stripChrPrefix(region.segment.chromosome) == "Y"){ continue; }
        totalBaf += region.bafCount;
        if(region.observedBaf > 0.57 && region.bafCount > 1){ highBaf += region.bafCount; }
    }
    const bool aneuploid = totalBaf > 0 && greaterOrEqual(static_cast<double>(highBaf) / totalBaf, 0.008);
    const bool diploidHighPurity = lowest.purity > 0.92 && lowest.ploidy > 1.8 && lowest.ploidy < 2.2;
    const bool highlyDiploid = greaterOrEqual(score.maxDiploidProportion, 0.97);
    const std::string method = (!aneuploid && (diploidHighPurity || highlyDiploid)) ? "NO_TUMOR" : "NORMAL";
    if(method == "NO_TUMOR"){
        FittedPurity converted = lowest;
        converted.purity = 1.0;
        converted.ploidy = 2.0;
        for(const auto &candidate : fits){
            if(std::abs(candidate.ploidy - 2.0) < 0.001 && std::abs(candidate.purity - 1.0) < 0.001){
                converted.score = candidate.score;
                converted.diploidProportion = candidate.diploidProportion;
                converted.somaticPenalty = candidate.somaticPenalty;
                break;
            }
        }
        return BestFit{converted, score, method};
    }
    return BestFit{lowest, score, method};
}

std::vector<ObservedRegion> fitObservedRegions(
        const std::vector<ObservedRegion> &regions, const FittedPurity &fit, int averageTumorDepth){
    std::vector<ObservedRegion> result;
    result.reserve(regions.size());
    const double ambiguous = expectedBaf(averageTumorDepth);
    for(const auto &source : regions){
        if(lp::stripChrPrefix(source.segment.chromosome) == "Y"){ continue; }
        ObservedRegion region = source;
        region.tumorCopyNumber = adjustedCopyNumber(region.observedTumorRatio, 1, fit.purity, fit.normFactor);
        region.tumorBaf = impliedBaf(region.tumorCopyNumber, region.observedBaf, fit.purity, fit.normFactor, ambiguous);
        region.refNormalisedCopyNumber = adjustedCopyNumber(
                region.observedTumorRatio, region.observedNormalRatio, fit.purity, fit.normFactor);
        const double major = region.tumorBaf * region.tumorCopyNumber;
        const double minor = region.tumorCopyNumber - major;
        region.majorAlleleCopyNumberDeviation = majorDeviation(fit.purity, fit.normFactor, major);
        region.minorAlleleCopyNumberDeviation = minorDeviation(fit.purity, fit.normFactor, minor);
        region.eventPenalty = 1 + 0.4 * std::min(std::abs(major - 1) + std::abs(minor - 1),
                                                1 + std::abs(major - 2) + std::abs(minor - 2));
        region.deviationPenalty = (region.minorAlleleCopyNumberDeviation + region.majorAlleleCopyNumberDeviation) * region.observedBaf;
        result.push_back(std::move(region));
    }
    return result;
}

void dumpFittingRegions(const std::vector<const ObservedRegion *> &regions){
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(regions.size());
    for(const auto *x : regions){ const auto &s = x->segment; rows.push_back({s.chromosome, s.start, s.chromosome + "\t" + std::to_string(s.start) + "\t" + std::to_string(s.end) + "\t" + germlineStatusName(x->germlineStatus) + "\t" + std::to_string(x->bafCount) + "\t" + lp::CpDump::num(x->observedBaf) + "\t" + lp::CpDump::num(x->observedNormalRatio) + "\t" + lp::CpDump::num(x->observedTumorRatio) + "\t0"}); }
    lp::CpDump::write("CP-P4-fitting-regions", "chromosome\tstart\tend\tgermlineStatus\tbafCount\tobservedBAF\tobservedNormalRatio\tobservedTumorRatio\tvariantCount", rows);
}

void dumpPurityGrid(const std::vector<FittedPurity> &fits){
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(fits.size());
    for(const auto &x : fits){ rows.push_back({"", 0, lp::CpDump::num(x.purity) + "\t" + lp::CpDump::num(x.normFactor) + "\t" + lp::CpDump::num(x.ploidy) + "\t" + lp::CpDump::num(x.score) + "\t" + lp::CpDump::num(x.diploidProportion) + "\t0"}); }
    lp::CpDump::write("CP-P5-purity-grid", "purity\tnormFactor\tploidy\tscore\tdiploidProportion\tsomaticPenalty", rows);
}

void dumpBestFit(const BestFit &best){
    const auto &x = best.fit; const auto &s = best.score;
    std::vector<lp::CpDump::Row> rows{{"", 0, best.method + "\t" + lp::CpDump::num(x.purity) + "\t" +
        lp::CpDump::num(x.normFactor) + "\t" + lp::CpDump::num(x.ploidy) + "\t" + lp::CpDump::num(x.score) + "\t" +
        lp::CpDump::num(x.diploidProportion) + "\t" + lp::CpDump::num(x.somaticPenalty) + "\t" +
        lp::CpDump::num(s.minPurity) + "\t" + lp::CpDump::num(s.maxPurity) + "\t" + lp::CpDump::num(s.minPloidy) + "\t" +
        lp::CpDump::num(s.maxPloidy) + "\t" + lp::CpDump::num(s.minDiploidProportion) + "\t" + lp::CpDump::num(s.maxDiploidProportion)}};
    lp::CpDump::write("CP-P6-best-fit", "method\tpurity\tnormFactor\tploidy\tscore\tdiploidProportion\tsomaticPenalty\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion", rows);
}

void dumpFittedRegions(const std::vector<ObservedRegion> &regions){
    std::vector<lp::CpDump::Row> rows;
    rows.reserve(regions.size());
    for(const auto &x : regions){ const auto &s=x.segment; rows.push_back({s.chromosome,s.start,
        s.chromosome+"\t"+std::to_string(s.start)+"\t"+std::to_string(s.end)+"\t"+(s.ratioSupport?"true":"false")+"\t"+
        segmentSupportName(s.support)+"\t"+std::to_string(x.bafCount)+"\t"+lp::CpDump::num(x.observedBaf)+"\t"+
        std::to_string(x.depthWindowCount)+"\t"+lp::CpDump::num(x.observedTumorRatio)+"\t"+lp::CpDump::num(x.observedNormalRatio)+"\t"+
        lp::CpDump::num(x.unnormalisedObservedNormalRatio)+"\t"+germlineStatusName(x.germlineStatus)+"\t"+(s.svCluster?"true":"false")+"\t"+
        lp::CpDump::num(x.gcContent)+"\t"+std::to_string(s.minStart)+"\t"+std::to_string(s.maxStart)+"\t"+
        lp::CpDump::num(x.minorAlleleCopyNumberDeviation)+"\t"+lp::CpDump::num(x.majorAlleleCopyNumberDeviation)+"\t"+
        lp::CpDump::num(x.deviationPenalty)+"\t"+lp::CpDump::num(x.eventPenalty)+"\t"+lp::CpDump::num(x.refNormalisedCopyNumber)+"\t"+
        lp::CpDump::num(x.tumorCopyNumber)+"\t"+lp::CpDump::num(x.tumorBaf)+"\t"+lp::CpDump::num(x.fittedTumorCopyNumber)+"\t"+
        lp::CpDump::num(x.fittedBaf)}); }
    lp::CpDump::write("CP-P7-fitted-regions", "chromosome\tstart\tend\tratioSupport\tsupport\tbafCount\tobservedBAF\tdepthWindowCount\tobservedTumorRatio\tobservedNormalRatio\tunnormalisedObservedNormalRatio\tgermlineStatus\tsvCluster\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumberDeviation\tmajorAlleleCopyNumberDeviation\tdeviationPenalty\teventPenalty\trefNormalisedCopyNumber\ttumorCopyNumber\ttumorBAF\tfittedTumorCopyNumber\tfittedBAF", rows);
}

}
