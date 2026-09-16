#include "PurpleWriters.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>

#include "../common/HumanChromosome.h"

namespace purple {
namespace {

std::string joinPath(const std::string &directory, const std::string &name){
    return directory + (directory.empty() || directory.back() == '/' ? "" : "/") + name;
}

std::string number(double value){
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);
    return buffer;
}

std::ofstream output(const std::string &path){
    std::ofstream stream(path);
    if(!stream){ throw std::runtime_error("unable to open output: " + path); }
    return stream;
}

double minor(const PurpleCopyNumber &copyNumber){
    if(copyNumber.averageActualBaf < 0.5 - 1e-10){ return 0; }
    return std::max(0.0, (1 - copyNumber.averageActualBaf) * copyNumber.averageTumorCopyNumber);
}

double major(const PurpleCopyNumber &copyNumber){ return copyNumber.averageTumorCopyNumber - minor(copyNumber); }

void writePurity(const std::string &path, const InputData &inputs, const BestFit &bestFit, const SummaryContext &summary){
    auto out = output(path);
    out << "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\trunMode\ttargeted\n";
    out << number(bestFit.fit.purity) << '\t' << number(bestFit.fit.normFactor) << '\t' << number(bestFit.fit.score) << '\t'
        << number(bestFit.fit.diploidProportion) << '\t' << number(bestFit.fit.ploidy) << '\t' << genderName(inputs.amberGender)
        << '\t' << bestFit.method << '\t' << number(summary.polyclonalProportion) << '\t' << number(bestFit.score.minPurity)
        << '\t' << number(bestFit.score.maxPurity) << '\t' << number(bestFit.score.minPloidy) << '\t' << number(bestFit.score.maxPloidy)
        << '\t' << number(bestFit.score.minDiploidProportion) << '\t' << number(bestFit.score.maxDiploidProportion) << '\t'
        << number(bestFit.fit.somaticPenalty) << '\t' << (summary.wholeGenomeDuplication ? "true" : "false")
        << "\t0.0000\tUNKNOWN\t0\tUNKNOWN\t0.0000\tUNKNOWN\t0\tTUMOR\tfalse\n";
}

void writeRange(const std::string &path, const std::vector<FittedPurity> &fits){
    auto out = output(path);
    out << "purity\tnormFactor\tscore\tdiploidProportion\tploidy\tsomaticPenalty\n";
    for(const auto &fit : fits){
        out << number(fit.purity) << '\t' << number(fit.normFactor) << '\t' << number(fit.score) << '\t'
            << number(fit.diploidProportion) << '\t' << number(fit.ploidy) << '\t' << number(fit.somaticPenalty) << '\n';
    }
}

void writeSegments(const std::string &path, const std::vector<ObservedRegion> &regions){
    auto out = output(path);
    out << "chromosome\tstart\tend\tgermlineStatus\tsvCluster\tbafCount\tobservedBAF\tminorAlleleCopyNumberDeviation\tobservedTumorRatio\tobservedNormalRatio\tunnormalisedObservedNormalRatio\tmajorAlleleCopyNumberDeviation\tdeviationPenalty\ttumorCopyNumber\tfittedTumorCopyNumber\tfittedBAF\trefNormalisedCopyNumber\tratioSupport\tsupport\tdepthWindowCount\ttumorBAF\tgcContent\teventPenalty\tminStart\tmaxStart\n";
    for(const auto &region : regions){
        const auto &segment = region.segment;
        out << segment.chromosome << '\t' << segment.start << '\t' << segment.end << '\t' << germlineStatusName(region.germlineStatus)
            << '\t' << (segment.svCluster ? "true" : "false") << '\t' << region.bafCount << '\t' << number(region.observedBaf)
            << '\t' << number(region.minorAlleleCopyNumberDeviation) << '\t' << number(region.observedTumorRatio)
            << '\t' << number(region.observedNormalRatio) << '\t' << number(region.unnormalisedObservedNormalRatio)
            << '\t' << number(region.majorAlleleCopyNumberDeviation) << '\t' << number(region.deviationPenalty)
            << '\t' << number(region.tumorCopyNumber) << '\t' << number(region.fittedTumorCopyNumber) << '\t' << number(region.fittedBaf)
            << '\t' << number(region.refNormalisedCopyNumber) << '\t' << (segment.ratioSupport ? "true" : "false")
            << '\t' << segmentSupportName(segment.support) << '\t' << region.depthWindowCount << '\t' << number(region.tumorBaf)
            << '\t' << number(region.gcContent) << '\t' << number(region.eventPenalty) << '\t' << segment.minStart << '\t' << segment.maxStart << '\n';
    }
}

void writeCopyNumbers(const std::string &path, const std::vector<PurpleCopyNumber> &copyNumbers){
    auto out = output(path);
    out << "chromosome\tstart\tend\tcopyNumber\tbafCount\tobservedBAF\tbaf\tsegmentStartSupport\tsegmentEndSupport\tmethod\tdepthWindowCount\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumber\tmajorAlleleCopyNumber\n";
    for(const auto &copyNumber : copyNumbers){
        out << copyNumber.chromosome << '\t' << copyNumber.start << '\t' << copyNumber.end << '\t' << number(copyNumber.averageTumorCopyNumber)
            << '\t' << copyNumber.bafCount << '\t' << number(copyNumber.averageObservedBaf) << '\t' << number(copyNumber.averageActualBaf)
            << '\t' << segmentSupportName(copyNumber.segmentStartSupport) << '\t' << segmentSupportName(copyNumber.segmentEndSupport)
            << '\t' << copyNumberMethodName(copyNumber.method) << '\t' << copyNumber.depthWindowCount << '\t' << number(copyNumber.gcContent)
            << '\t' << copyNumber.minStart << '\t' << copyNumber.maxStart << '\t' << number(minor(copyNumber)) << '\t' << number(major(copyNumber)) << '\n';
    }
}

const std::array<int,23> CENTROMERE_38 = {123605523,93139351,92214016,50726026,48272854,59191911,59498944,44955505,44377363,40640102,52751711,35977330,17025624,17086762,18362627,37295920,24849830,18161053,25844927,28237290,11890184,14004553,60509061};

void emitArm(std::ofstream &out, const std::string &chromosome, char arm, std::vector<const PurpleCopyNumber *> values){
    if(values.empty()){ return; }
    long long totalLength = 0;
    double weighted = 0;
    double minimum = std::numeric_limits<double>::max();
    // Java Double.MIN_VALUE is the smallest positive subnormal, not the most
    // negative value. Preserve that behavior for all-negative arms.
    double maximum = std::numeric_limits<double>::denorm_min();
    for(const auto *value : values){
        const int length = value->end - value->start + 1;
        totalLength += length;
        weighted += static_cast<double>(length) * value->averageTumorCopyNumber;
        minimum = std::min(minimum, value->averageTumorCopyNumber);
        maximum = std::max(maximum, value->averageTumorCopyNumber);
    }
    std::sort(values.begin(), values.end(), [](const auto *left, const auto *right){
        if(left->averageTumorCopyNumber != right->averageTumorCopyNumber){ return left->averageTumorCopyNumber < right->averageTumorCopyNumber; }
        return left->start < right->start;
    });
    const int halfLength = static_cast<int>(totalLength / 2);
    int runningLength = 0;
    double median = 0;
    for(const auto *value : values){
        runningLength += value->end - value->start + 1;
        if(runningLength > halfLength){ median = value->averageTumorCopyNumber; break; }
    }
    out << chromosome << '\t' << arm << '\t' << number(weighted / totalLength) << '\t' << number(median)
        << '\t' << number(minimum) << '\t' << number(maximum) << '\n';
}

void writeArms(const std::string &path, const std::vector<PurpleCopyNumber> &copyNumbers){
    auto out = output(path);
    out << "chromosome\tarm\tmeanCopyNumber\tmedianCopyNumber\tminCopyNumber\tmaxCopyNumber\n";
    for(int chromosome = 1; chromosome <= 23; ++chromosome){
        const std::string name = chromosome == 23 ? "X" : std::to_string(chromosome);
        std::vector<const PurpleCopyNumber *> p, q;
        for(const auto &copyNumber : copyNumbers){
            if(lp::stripChrPrefix(copyNumber.chromosome) != name){ continue; }
            (copyNumber.start < CENTROMERE_38[chromosome - 1] ? p : q).push_back(&copyNumber);
        }
        const bool acrocentric = chromosome == 13 || chromosome == 14 || chromosome == 15 || chromosome == 21 || chromosome == 22;
        if(!acrocentric){ emitArm(out, name, 'P', std::move(p)); }
        emitArm(out, name, 'Q', std::move(q));
    }
}

void writeQc(const std::string &path, const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary){
    auto out = output(path);
    out << "QCStatus\t" << summary.qcStatus << "\nMethod\t" << bestFit.method
        << "\nCopyNumberSegments\t" << copyNumbers.size() << "\nUnsupportedCopyNumberSegments\t0\nPurity\t" << number(bestFit.fit.purity)
        << "\nAmberGender\t" << genderName(inputs.amberGender) << "\nCobaltGender\t" << genderName(inputs.cobaltGender)
        << "\nDeletedGenes\t" << summary.deletedGenes << "\nContamination\t" << number(inputs.contamination)
        << "\nGermlineAberrations\tNONE\nAmberMeanDepth\t" << inputs.averageTumorDepth << "\nLohPercent\t" << number(summary.lohPercent)
        << "\nTincLevel\t0.0000\n";
}

}

void writeCoreOutputs(const std::string &outputDir, const InputData &inputs,
        const std::vector<FittedPurity> &fits, const BestFit &bestFit,
        const std::vector<ObservedRegion> &fittedRegions,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary){
    if(outputDir.empty()){ return; }
    mkdir(outputDir.c_str(), 0755);
    const std::string prefix = inputs.sampleId + ".purple.";
    writePurity(joinPath(outputDir, prefix + "purity.tsv"), inputs, bestFit, summary);
    writeRange(joinPath(outputDir, prefix + "purity.range.tsv"), fits);
    writeSegments(joinPath(outputDir, prefix + "segment.tsv"), fittedRegions);
    writeCopyNumbers(joinPath(outputDir, prefix + "cnv.somatic.tsv"), copyNumbers);
    writeArms(joinPath(outputDir, prefix + "chromosome_arm.tsv"), copyNumbers);
    writeQc(joinPath(outputDir, prefix + "qc"), inputs, bestFit, copyNumbers, summary);
}

}
