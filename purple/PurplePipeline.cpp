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

namespace purple {

InputData loadInputs(const PipelineConfig &cfg)
{
    return loadTumorOnlyInputs(cfg.sampleId, cfg.amberDir, cfg.cobaltDir);
}

void runFromInputs(const PipelineConfig &cfg, const InputData &inputs,
        const ChromosomeLengths &lengths, const LohFillInput *lohFill)
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
    const auto fittingRegions = selectFittingRegions(observed);
    dumpFittingRegions(fittingRegions);
    const auto fits = fitPurityGrid(fittingRegions, inputs.averageTumorDepth, cfg.threads);
    dumpPurityGrid(fits);
    const auto bestFit = selectTumorOnlyBestFit(fits, observed);
    dumpBestFit(bestFit);
    const auto fittedRegions = fitObservedRegions(observed, bestFit.fit, inputs.averageTumorDepth,
            inputs.cobaltGender);
    dumpFittedRegions(fittedRegions);
    const auto copyNumbers = buildCopyNumbers(fittedRegions, bestFit.fit, inputs.cobaltGender);
    dumpCopyNumbers(copyNumbers);
    const auto summary = buildSummaryContext(inputs, bestFit, copyNumbers, cfg.ensemblDataDir);
    dumpSummaryContext(inputs, bestFit, copyNumbers, summary);
    writeCoreOutputs(cfg.outputDir, inputs, fits, bestFit, fittedRegions, copyNumbers, summary);
    // 放在 writeCoreOutputs 之後：輸出目錄由它建立。
    if(lohFill != nullptr){ writeLohFillTsv(cfg.outputDir, inputs.sampleId, fill); }
    std::cerr << "PURPLE P10 core-output stage complete: " << inputs.bafs.size() << " BAF, "
              << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
              << inputs.cobaltTumorPcf.size() << " Cobalt PCF, " << segments.size() << " support segments, "
              << observed.size() << " observed regions, " << fits.size() << " purity/ploidy candidates\n";
}

void run(const PipelineConfig &cfg, const LohFillInput *lohFill)
{
    runFromInputs(cfg, loadInputs(cfg), readChromosomeLengths(cfg.refGenome), lohFill);
}

}
