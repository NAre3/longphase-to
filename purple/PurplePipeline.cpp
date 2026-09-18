#include "PurplePipeline.h"

#include <algorithm>
#include <iostream>

#include "PurpleInput.h"
#include "PurpleSegmentation.h"
#include "PurpleObserved.h"
#include "PurpleFitting.h"
#include "PurpleCopyNumber.h"
#include "PurpleSummary.h"
#include "PurpleWriters.h"

namespace purple {

InputData loadInputs(const PipelineConfig &cfg)
{
    return loadTumorOnlyInputs(cfg.sampleId, cfg.amberDir, cfg.cobaltDir);
}

void runFromInputs(const PipelineConfig &cfg, const InputData &inputs)
{
    // 這個函式的內容原本是 PurpleApplication.cpp 的 main() 主體（34b8d97）。
    // 抽出時只改了取值來源（argv -> cfg），沒有改動任何呼叫順序或引數。
    dumpInputCheckpoint(inputs);
    const auto segments = createSupportSegments(inputs, cfg.refGenome);
    dumpSupportSegments(segments);
    const auto observed = createObservedRegions(inputs, segments);
    dumpObservedRegions(observed);
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
    std::cerr << "PURPLE P10 core-output stage complete: " << inputs.bafs.size() << " BAF, "
              << inputs.ratios.size() << " ratio, " << inputs.amberPcf.size() << " Amber PCF, "
              << inputs.cobaltTumorPcf.size() << " Cobalt PCF, " << segments.size() << " support segments, "
              << observed.size() << " observed regions, " << fits.size() << " purity/ploidy candidates\n";
}

void run(const PipelineConfig &cfg)
{
    runFromInputs(cfg, loadInputs(cfg));
}

}
