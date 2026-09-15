#include "PhasingProcess.h"
#include "PhasingGraph.h"
#include "ParsingBam.h"
#include "SomaticRefinementPolicy.h"

#include <set>

#include "amber/AmberPipeline.h"
#include "amber/SharedScanSink.h"
#include "common/CpDump.h"

PhasingProcess::PhasingProcess(PhasingParameters params)
{
    std::cerr<< "LongPhase-TO Ver " << params.version << "\n";
    std::cerr<< "\n";
    std::cerr<< "--- File Parameter --- \n";
    std::cerr<< "SNP File       : " << params.snpFile       << "\n";
    std::cerr<< "SV  File       : " << params.svFile        << "\n";
    std::cerr<< "MOD File       : " << params.modFile       << "\n";
    std::cerr<< "PON File       : " << params.ponFile       << "\n";
    std::cerr<< "Strict PON File: " << params.strictPonFile << "\n";
    std::cerr<< "REF File       : " << params.fastaFile     << "\n";
    std::cerr<< "Output Prefix  : " << params.resultPrefix  << "\n";
    std::cerr<< "Generate Dot   : " << ( params.generateDot ? "True" : "False" ) << "\n";
    std::cerr<< "Output LOH     : " << ( params.outputLOH ? "True" : "False" ) << "\n";
    std::cerr<< "Output SGE     : " << ( params.outputSGE ? "True" : "False" ) << "\n";
    std::cerr<< "Output LGE     : " << ( params.outputLGE ? "True" : "False" ) << "\n";
    std::cerr<< "Output GE      : " << ( params.outputGE ? "True" : "False" ) << "\n";
    std::cerr<< "BAM File       : ";
    for( auto file : params.bamFile){
        std::cerr<< file <<" " ;
    }
    std::cerr << "\n";

    std::cerr<< "\n";
    std::cerr<< "--- Phasing Parameter --- \n";
    std::cerr<< "Caller             : " << params.callerStr << "\n";
    std::cerr<< "Phase Indel        : " << ( params.phaseIndel ? "True" : "False" )  << "\n";
    std::cerr<< "Distance Threshold : " << params.distance        << "\n";
    std::cerr<< "Connect Adjacent   : " << params.connectAdjacent << "\n";
    std::cerr<< "Somatic Connect Adjacent   : " << params.somaticConnectAdjacent << "\n";
    std::cerr<< "Edge Threshold     : " << params.edgeThreshold   << "\n";
    std::cerr<< "Overlap Threshold  : " << params.overlapThreshold   << "\n";
    std::cerr<< "Mapping Quality    : " << params.mappingQuality  << "\n";
    std::cerr<< "Mismatch Rate      : " << params.mismatchRate  << "\n";
    std::cerr<< "Variant Confidence : " << params.snpConfidence   << "\n";
    std::cerr<< "ReadTag Confidence : " << params.readConfidence  << "\n";
    std::cerr<< "\n";

    std::cerr<< "--- Methylation XGBoost Parameter --- \n";
    std::cerr<< "Methyl XGBoost    : " << ( params.enableMethylXgb ? "True" : "False" ) << "\n";
    if(params.enableMethylXgb){
        std::cerr<< "Model Source      : embedded C++\n";
        std::cerr<< "SNV/indel Cutoff : " << params.methylXgbSnvThreshold
                 << "/" << params.methylXgbIndelThreshold << "\n";
        std::cerr<< "Purity Eligibility: <= " << somatic_refinement::kMethylXgbMaxPurity << "\n";
        std::cerr<< "Methyl Window     : " << params.methylXgbWindow << "\n";
        std::cerr<< "Meth High/Low     : " << params.methylXgbMethHigh << "/" << params.methylXgbMethLow << "\n";
    }
    std::cerr<< "\n";

    std::time_t processBegin = time(NULL);

    // load SNP vcf file
    std::time_t begin = time(NULL);
    std::cerr<< "parsing VCF ... ";
    SnpParser snpFile(params);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load PON vcf and strict PON vcf file
    begin = time(NULL);
    std::cerr<< "parsing PON VCF ... ";
    snpFile.setGermline(params.ponFile, params.strictPonFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // load SV vcf file
    begin = time(NULL);
    std::cerr<< "parsing SV VCF ... ";
    SVParser svFile(params, snpFile);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    //Parse mod vcf file
	begin = time(NULL);
	std::cerr<< "parsing Meth VCF ... ";
    METHParser modFile(params, snpFile, svFile);
	std::cerr<< difftime(time(NULL), begin) << "s\n";

    // parsing ref fasta
    begin = time(NULL);
    std::cerr<< "reading reference ... ";
    std::vector<int> last_pos;
    for(auto chr :snpFile.getChrVec()){
        last_pos.push_back(snpFile.getLastSNP(chr));
    }
    FastaParser fastaParser(params.fastaFile, snpFile.getChrVec(), last_pos, params.numThreads);
    std::cerr<< difftime(time(NULL), begin) << "s\n";

    // get all detected chromosome
    std::vector<std::string> chrName = snpFile.getChrVec();

    // record chromosome info
    std::map<std::string, ChrInfo> chrInfoMap;
    // Initialize chrInfoMap entries for each chromosome name in advance
    // to avoid potential issues during parallel processing later
    for (std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end(); chrIter++)    {
        chrInfoMap[*chrIter] = ChrInfo();
    }

    // init data structure and get core n
    htsThreadPool threadPool = {NULL, 0};

    // creat thread pool
    if (!(threadPool.pool = hts_tpool_init(params.numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    if (!params.disableCalling){
        std::cerr << "parsing BAM, somatic calling, and phasing" << std::endl;
    }else{
        std::cerr << "parsing BAM, and phasing" << std::endl;
    }
    begin = time(NULL);

    // ---- AMBER prescan（EXP-I02）----
    //
    // CP-A1 → CP-A2b 全部在共用走訪之前完成，產出 evidence 與 RegionTask。
    // 之後每條染色體各自建一個 ContigSink，由共用走訪把 read 路由進去。
    amber::PipelineConfig amberCfg;
    amber::PrescanResult amberPrescan;
    std::vector<amber::ContigSink *> amberSinks(chrName.size(), nullptr);
    std::vector<std::string> scanContigs = chrName;

    if(params.amberEnabled()){
        amberCfg.lociPath = params.amberLoci;
        amberCfg.bedPath = params.amberExcludedBed;
        amberCfg.tumorBam = params.bamFile.front();
        amberCfg.outputDir = params.amberOutputDir;
        amberCfg.sampleId = params.amberSampleId;
        amberCfg.minBaseQuality = params.amberMinBaseQuality;
        amberCfg.minMappingQuality = params.amberMinMapQuality;
        amberCfg.threads = params.numThreads;

        lp::CpDump::setDir(params.amberCpDumpDir);

        begin = time(NULL);
        std::cerr<< "AMBER prescan ... ";
        amberPrescan = amber::prescan(amberCfg);
        std::cerr<< difftime(time(NULL), begin) << "s\n";

        // 掃描的 contig 清單 = 候選 VCF 的 contig ∪ AMBER 的 contig（C4／D7）。
        // 候選 VCF 沒有的 contig 仍必須被掃描，否則該 contig 上的 AMBER 位點整條缺席。
        const std::vector<std::pair<std::string, std::pair<std::size_t, std::size_t>>> taskIndex =
                amber::indexTasksByChromosome(amberPrescan.tasks);

        std::set<std::string> vcfContigs(chrName.begin(), chrName.end());
        for(const auto &entry : taskIndex){
            if(vcfContigs.find(entry.first) == vcfContigs.end()){
                scanContigs.push_back(entry.first);
            }
        }

        // sink 的擁有權放在這個 vector，索引與 scanContigs 對齊。
        // 先一次配置完，迴圈內不再有任何容器改動 → 平行迴圈無共享寫入。
        amberSinks.assign(scanContigs.size(), nullptr);
        for(const auto &entry : taskIndex){
            const std::size_t slot = std::find(scanContigs.begin(), scanContigs.end(), entry.first)
                    - scanContigs.begin();
            if(slot >= scanContigs.size()){
                continue;
            }
            amberSinks[slot] = new amber::ContigSink(
                    amberPrescan.tasks.data() + entry.second.first,
                    amberPrescan.tasks.data() + entry.second.second,
                    amberCfg.minMappingQuality, amberCfg.minBaseQuality);
        }
    }

    // loop all chromosome
    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(int scanIdx = 0; scanIdx < static_cast<int>(scanContigs.size()); scanIdx++ ){

        // 取副本：VairiantGraph 的建構子要的是 std::string&（非 const）
        std::string contigName = scanContigs[scanIdx];
        std::time_t chrbegin = time(NULL);
        amber::ContigSink *amberSink = amberSinks[scanIdx];

        // get lase SNP variant position
        // scanContigs 的尾段可能是 AMBER 專屬的 contig，那些不在 snpFile 內。
        int lastSNPpos = (scanIdx < static_cast<int>(chrName.size()))
                ? snpFile.getLastSNP(contigName) : -1;

        // 共用走訪的右界（design.md D9）
        int scanRightEdge = lastSNPpos;
        if( amberSink != nullptr ){
            scanRightEdge = std::max(scanRightEdge, amberSink->maxTaskEnd());
        }

        // therer is no variant on SNP file.
        if( lastSNPpos == -1 ){
            // D7：沒有候選變異的 contig 現在仍必須被掃描給 AMBER，只跳過 clip／graph／phasing。
            // BamParser 在此建不出來（D4：建構子對空 variant map 直接 exit(1)），
            // 故走另一個進入點。
            if( amberSink != nullptr && scanRightEdge > 0 ){
                amberOnlyContigScan(params.bamFile.front(), contigName, scanRightEdge,
                        threadPool, params, *amberSink);
            }
            continue;
        }

        ChrInfo &chrInfo = chrInfoMap[contigName];

	    // fetch chromosome string
        std::string chr_reference = fastaParser.chrString.at(contigName);
        // create a bam parser object and prepare to fetch varint from each vcf file
	    BamParser *bamParser = new BamParser(contigName, params.bamFile, snpFile, svFile, modFile, chr_reference);
        // use to store variant
        std::vector<ReadVariant> *readVariantVec = new std::vector<ReadVariant>();
        // run fetch variant process
        bamParser->direct_detect_alleles(lastSNPpos, scanRightEdge, threadPool, params, *readVariantVec, chrInfo.clipCount, chr_reference, amberSink);
        // free memory
        delete bamParser;

        const bool hasPhasingVariant = std::any_of(
            readVariantVec->begin(),
            readVariantVec->end(),
            [](const ReadVariant &read) { return !read.variantVec.empty(); }
        );
        // BAM files are partial or no read supports a phasing variant on this chromosome.
        if(!hasPhasingVariant){
            delete readVariantVec;
            continue;
        }

        // create a clip object and prepare to detect Interval
        Clip *clip = new Clip(contigName);
        // get the interval of the genomic event
        clip->detectGenomicEventInterval(chrInfo.clipCount, chrInfo.largeGenomicEventInterval, chrInfo.smallGenomicEventRegion);
        // get the region of the LOH
        clip->detectLOHRegion(snpFile, chrInfo.LOHSegments);
        // free memory
        delete clip;

        // create a graph object and prepare to phasing.
        VairiantGraph *vGraph = new VairiantGraph(chr_reference, params, contigName);
        chrInfo.vGraph = vGraph;
        // trans read-snp info to edge info
        vGraph->addEdge(readVariantVec, chrInfo.LOHSegments);
        if(!params.disableCalling){
            // run somatic calling algorithm
            vGraph->somaticCalling(snpFile.getVariants(contigName));
        }else{
            vGraph->tagSomatic(snpFile.getVariants(contigName));
        }
        // run main algorithm
        vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, &chrInfo.ploidyRatioMap);
        std::cerr<< "(" << contigName << "," << difftime(time(NULL), chrbegin) << "s)";
    }

    // ---- AMBER postscan（EXP-I02）----
    //
    // CP-A3 → CP-A9 與三個 stage 輸出。呼叫的是 amber_port 用的同一組函式。
    if(params.amberEnabled()){
        std::cerr << std::endl;

        // 每條染色體一個 sink，各自持有自己的計數。兩個計數器都是整數加總，
        // 與相加順序無關（BamEvidenceReader.cpp:284-289 的同一個論證）。
        amber::BamScanStats amberStats;
        for(amber::ContigSink *sink : amberSinks){
            if(sink == nullptr){
                continue;
            }
            amberStats.recordsConsumed += sink->stats().recordsConsumed;
            amberStats.nonAcgtnBases += sink->stats().nonAcgtnBases;
        }

        // recordsConsumed 的語義與 amber_port 不同（design.md D3）：此處是「通過 slicer
        // filter 的 read 數」，amber_port 是「per-region 造訪次數」。刻意不重現，
        // 該欄位不進任何 checkpoint 或 stage 輸出。
        std::fprintf(stderr, "AMBER shared scan: consumed %llu reads (per-read, not per-region), "
                "non-ACGTN bases at evaluated positions: %llu\n",
                static_cast<unsigned long long>(amberStats.recordsConsumed),
                static_cast<unsigned long long>(amberStats.nonAcgtnBases));

        for(amber::ContigSink *sink : amberSinks){
            delete sink;
        }
        amberSinks.clear();

        begin = time(NULL);
        std::cerr<< "AMBER postscan ... ";
        amber::postscan(amberCfg, amberPrescan.evidence, amberStats);
        std::cerr<< difftime(time(NULL), begin) << "s\n";
    }

    std::map<std::string, std::map<double, int>> mergedPloidyRatioMap;
    for (const auto& chrInfo : chrInfoMap) {
        mergedPloidyRatioMap[chrInfo.first] = chrInfo.second.ploidyRatioMap;
    }
    double purity = params.purity;
    if (purity < 0.0) {
        purity = PurityCalculator::getPurity(mergedPloidyRatioMap, params.resultPrefix, params.caller, chrInfoMap, fastaParser.chrLength);
    }
    std::cerr << std::endl;
    std::cerr << "purity: " << purity << std::endl;
    const somatic_refinement::Plan refinementPlan =
        somatic_refinement::makeResolvedPlan(purity, params.enableMethylXgb);
    const bool highPurity = refinementPlan.convertNonGermlineToSomatic;
    const bool runMethylXgb = refinementPlan.runMethylXgb;
    if(runMethylXgb){
        std::cerr << "run embedded methyl XGBoost refinement ... ";
        std::vector<MethylXgbApplySummary> methylXgbSummaries(chrName.size());
        #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
        for(int chrIdx = 0; chrIdx < static_cast<int>(chrName.size()); chrIdx++){
            ChrInfo &chrInfo = chrInfoMap[chrName[chrIdx]];
            if(chrInfo.vGraph == nullptr) {
                continue;
            }
            methylXgbSummaries[chrIdx] = chrInfo.vGraph->refineSomaticWithMethylXgb();
        }

        MethylXgbApplySummary methylXgbSummary;
        for(const auto &summary : methylXgbSummaries){
            methylXgbSummary.applied += summary.applied;
            methylXgbSummary.somatic += summary.somatic;
            methylXgbSummary.nonSomatic += summary.nonSomatic;
        }
        std::cerr << methylXgbSummary.applied << " predictions, "
                  << methylXgbSummary.somatic << " somatic, "
                  << methylXgbSummary.nonSomatic << " non-somatic\n";
    }
    else if(params.enableMethylXgb){
        std::cerr << "skip methyl XGBoost refinement: purity " << purity
                  << " > " << somatic_refinement::kMethylXgbMaxPurity << "\n";
    }
    if(highPurity){
        std::cerr << "convert non-PON variants to somatic for high-purity sample" << std::endl;
    }
    if(refinementPlan.rerunPhasing){
        std::cerr << "second round phasing, ";
    }
    std::cerr << "export phasing result" << std::endl;

    #pragma omp parallel for schedule(dynamic) num_threads(params.numThreads)
    for(std::vector<std::string>::iterator chrIter = chrName.begin(); chrIter != chrName.end() ; chrIter++ ){
        std::time_t chrbegin = time(NULL);
        ChrInfo &chrInfo = chrInfoMap[*chrIter];
        if(chrInfo.vGraph == nullptr)continue;

        VairiantGraph *vGraph = chrInfo.vGraph;
        if(runMethylXgb){
            // reset phasing result
            chrInfo.posPhasingResult = PosPhasingResult();
            // run main algorithm after methyl XGBoost somatic refinement
            vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, nullptr);
        }
        else if(highPurity){
            // convert non-pon variants to somatic variants
            vGraph->convertNonGermlineToSomatic();
            // reset phasing result
            chrInfo.posPhasingResult = PosPhasingResult();
            // run main algorithm
            vGraph->phasingProcess(chrInfo.posPhasingResult, chrInfo.LOHSegments, nullptr);
        }
        // export phasing result
        vGraph->exportPhasingResult(chrInfo.posPhasingResult, chrInfo.LOHSegments);
        // generate dot file
        if(params.generateDot){
            vGraph->writingDotFile((*chrIter));
        }

        // release the memory used by the object.
        vGraph->destroy();
        delete vGraph;

        std::cerr<< "(" << (*chrIter) << "," << difftime(time(NULL), chrbegin) << "s)";
    }
    hts_tpool_destroy(threadPool.pool);

    // Transfer phasing results from chrInfoMap to chrPhasingResult
    ChrPhasingResult chrPhasingResult;
    for (const auto& chrInfo : chrInfoMap) {
        chrPhasingResult[chrInfo.first] = chrInfo.second.posPhasingResult;
    }

    std::cerr<< "\nparsing total:  " << difftime(time(NULL), begin) << "s\n";

    // write result to file
    GenomicWriter genomicWriter(params.resultPrefix, chrName, chrInfoMap);
    genomicWriter.measureTime("LOH", params.outputLOH, [&]() { genomicWriter.writeLOH(); });
    genomicWriter.measureTime("SGE", params.outputSGE, [&]() { genomicWriter.writeSGE(); });
    genomicWriter.measureTime("LGE", params.outputLGE, [&]() { genomicWriter.writeLGE(); });
    genomicWriter.measureTime("GE", params.outputGE, [&]() { genomicWriter.writeAllEvents(); });
    genomicWriter.measureTime("SNP", true, [&]() { snpFile.writeResult(chrPhasingResult, purity); });
    genomicWriter.measureTime("SV", params.svFile != "", [&]() { svFile.writeResult(chrPhasingResult); });
    genomicWriter.measureTime("MOD", params.modFile != "", [&]() { modFile.writeResult(chrPhasingResult); });

    std::cerr<< "\ntotal process: " << difftime(time(NULL), processBegin) << "s\n";

    return;
};

PhasingProcess::~PhasingProcess(){
};

