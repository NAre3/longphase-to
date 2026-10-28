#include "PurplePipeline.h"

#include <algorithm>
#include <iostream>

#include "PurpleInput.h"
#include "PurpleSegmentation.h"
#include "PurpleObserved.h"
#include "PurpleLohFill.h"
#include "PurpleFitting.h"
#include "PurpleCopyNumber.h"
#include "PurpleSummary.h"
#include "PurpleWriters.h"
#include "PurplePhasedBaf.h"
#include "PurpleLohTdf.h"
#include "PurpleScaleCheck.h"
#include "common/HumanChromosome.h"

#include <set>

namespace purple {

InputData loadInputs(const PipelineConfig &cfg)
{
    return loadTumorOnlyInputs(cfg.sampleId, cfg.amberDir, cfg.cobaltDir);
}

BestFit runFromInputs(const PipelineConfig &cfg, const InputData &inputs,
        const ChromosomeLengths &lengths, const LohFillInput *lohFill, const PhasedSnvMap *phased)
{
    // 這個函式的內容原本是 PurpleApplication.cpp 的 main() 主體。
    // 抽出時只改了取值來源（argv -> cfg），沒有改動任何呼叫順序或引數。
    // 之後在 observed regions 與擬合之間加入了 LOH 補回步驟（lohFill 非 nullptr 時）。
    dumpInputCheckpoint(inputs);
    const auto segments = createSupportSegments(inputs, lengths);
    dumpSupportSegments(segments);
    auto observed = createObservedRegions(inputs, segments);
    dumpObservedRegions(observed);
    // LOH 補回（PurpleLohFill.h）。CP-P3 是補前的值；之後的擬合、選解、copy number
    // 與輸出全部吃補後的值。lohFill 為 nullptr 時不做任何事，行為與加入補回之前相同。
    LohFillResult fill;
    if(lohFill != nullptr){
        fill = applyLohFill(observed, *lohFill);
        dumpLohFill(fill);
        std::cerr << "PURPLE LOH fill: retention " << fill.ret << " over " << fill.retSegments
                  << " non-LOH segments, filled " << fill.filled.size() << " segments\n";
    }
    auto fittingRegions = selectFittingRegions(observed, cfg.excludeAllosomes);
    dumpFittingRegions(fittingRegions);
    // phased 段 BAF：在 observed 的副本上改寫擬合段的 observedBaf。CP-P4 是改寫前的值。
    // aneuploid 判斷（selectTumorOnlyBestFit）仍用原 observed：它的 0.57 門檻是針對折疊後的 BAF 定的。
    std::vector<ObservedRegion> observedPhased;
    std::vector<char> halfAmbiguous;
    std::vector<ObservedRegion *> targets;
    PhasedBafResult phasedResult;
    if(phased != nullptr){
        observedPhased = observed;
        halfAmbiguous.assign(observed.size(), 0);
        for(const ObservedRegion *region : fittingRegions){
            const size_t index = static_cast<size_t>(region - observed.data());
            halfAmbiguous[index] = 1;
            targets.push_back(&observedPhased[index]);
        }
        std::set<std::pair<std::string, int>> filledKeys;
        for(const auto &record : fill.filled){ filledKeys.insert({lp::stripChrPrefix(record.chromosome), record.start}); }
        phasedResult = applyPhasedBaf(targets, inputs.bafs, *phased, filledKeys);
        dumpPhasedBaf(targets, phasedResult);
        std::cerr << "PURPLE phased BAF: kappa " << phasedResult.kappa << ", switch " << phasedResult.switchProbability
                  << ", " << phasedResult.regionsUpdated << " of " << targets.size() << " fitting regions updated\n";
        fittingRegions.assign(targets.begin(), targets.end());
    }
    const auto fits = fitPurityGrid(fittingRegions, inputs.averageTumorDepth, cfg.threads, phased != nullptr ? 0.5 : -1);
    dumpPurityGrid(fits);
    // LOH 段的 tumor DNA fraction 限制候選（PurpleLohTdf.h）。purity.range.tsv 仍寫完整格點；
    // 選解（含 NO_TUMOR 判斷與 purity.tsv 的分數範圍欄）只看範圍內的候選。
    LohTdfEstimate lohTdf;
    std::vector<FittedPurity> windowFits;
    const bool useLohTdf = phased != nullptr && cfg.lohTdfWindow;
    if(useLohTdf){
        lohTdf = estimateLohTdf(targets, phasedResult);
        if(lohTdf.valid){ windowFits = restrictToTdfWindow(fits, lohTdf.f0, cfg.lohTdfHalfWidth); }
        std::cerr << "PURPLE LOH TDF window: "
                  << (lohTdf.valid ? "f0 " + std::to_string(lohTdf.f0) : std::string("no segments"))
                  << " from " << lohTdf.segments.size() << " segments, " << windowFits.size() << " of "
                  << fits.size() << " candidates within " << cfg.lohTdfHalfWidth
                  << (windowFits.empty() ? " (not applied)" : "") << "\n";
    }
    auto bestFit = selectTumorOnlyBestFit(windowFits.empty() ? fits : windowFits, observed);
    dumpBestFit(bestFit);
    const auto fitRegions = [&](const FittedPurity &fit){
        return phased != nullptr
                ? fitObservedRegions(observedPhased, fit, inputs.averageTumorDepth, inputs.cobaltGender, &halfAmbiguous)
                : fitObservedRegions(observed, fit, inputs.averageTumorDepth, inputs.cobaltGender);
    };
    auto fittedRegions = fitRegions(bestFit.fit);
    dumpFittedRegions(fittedRegions);
    auto copyNumbers = buildCopyNumbers(fittedRegions, bestFit.fit, inputs.cobaltGender);
    dumpCopyNumbers(copyNumbers);
    // 段格點尺度檢查（PurpleScaleCheck.h）：判 ×2 且 ψ ≥ 2 時改選 2ψ 附近的格點候選，
    // 之後的 fitted regions、copy number 與所有輸出都用改後的解；CP-P6～P8 以改後的值覆寫。
    // 分數範圍欄（purity.tsv 的 min／max）沿用原選解的值。
    ScaleCheck scaleCheck;
    const FittedPurity originalFit = bestFit.fit;
    if(cfg.scaleCorrection){
        scaleCheck = checkScale(copyNumbers, bestFit.fit, fits, useLohTdf && lohTdf.valid, lohTdf.f0);
        std::cerr << "PURPLE scale check: s* " << scaleCheck.sMax << ", A " << scaleCheck.aMax << " over "
                  << scaleCheck.segments << " segments, " << (scaleCheck.computed ? scaleCheck.decision : std::string("NA"))
                  << (scaleCheck.applied ? ", ploidy " + std::to_string(bestFit.fit.ploidy) + " -> " + std::to_string(scaleCheck.candidate.ploidy)
                                         : std::string(", not applied")) << "\n";
        if(scaleCheck.applied){
            bestFit.fit = scaleCheck.candidate;
            dumpBestFit(bestFit);
            fittedRegions = fitRegions(bestFit.fit);
            dumpFittedRegions(fittedRegions);
            copyNumbers = buildCopyNumbers(fittedRegions, bestFit.fit, inputs.cobaltGender);
            dumpCopyNumbers(copyNumbers);
        }
    }
    const auto summary = buildSummaryContext(inputs, bestFit, copyNumbers, cfg.ensemblDataDir);
    dumpSummaryContext(inputs, bestFit, copyNumbers, summary);
    writeCoreOutputs(cfg.outputDir, inputs, fits, bestFit, fittedRegions, copyNumbers, summary);
    // 放在 writeCoreOutputs 之後：輸出目錄由它建立。
    if(lohFill != nullptr){ writeLohFillTsv(cfg.outputDir, inputs.sampleId, fill); }
    if(cfg.scaleCorrection){ writeScaleCheckTsv(cfg.outputDir, inputs.sampleId, scaleCheck, originalFit); }
    if(useLohTdf){ writeLohTdfTsv(cfg.outputDir, inputs.sampleId, lohTdf, cfg.lohTdfHalfWidth, windowFits.size(), fits.size()); }
    std::cerr << "PURPLE core outputs written: " << inputs.bafs.size() << " BAF, "
              << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
              << inputs.cobaltTumorPcf.size() << " Cobalt PCF, " << segments.size() << " support segments, "
              << observed.size() << " observed regions, " << fits.size() << " purity/ploidy candidates\n";
    return bestFit;
}

void run(const PipelineConfig &cfg, const LohFillInput *lohFill)
{
    runFromInputs(cfg, loadInputs(cfg), readChromosomeLengths(cfg.refGenome), lohFill);
}

}
