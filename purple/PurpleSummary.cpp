#include "PurpleSummary.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

constexpr double EPS = 1e-10;

std::string joinPath(const std::string &directory, const std::string &name){
    return directory + (directory.empty() || directory.back() == '/' ? "" : "/") + name;
}

std::vector<std::string> csvFields(const std::string &line){
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for(std::size_t i = 0; i < line.size(); ++i){
        const char c = line[i];
        if(c == '"'){
            if(quoted && i + 1 < line.size() && line[i + 1] == '"'){ field.push_back('"'); ++i; }
            else{ quoted = !quoted; }
        }else if(c == ',' && !quoted){ fields.push_back(field); field.clear(); }
        else{ field.push_back(c); }
    }
    fields.push_back(field);
    return fields;
}

double minorAlleleCopyNumber(const PurpleCopyNumber &copyNumber){
    if(copyNumber.averageActualBaf < 0.5 - EPS){ return 0; }
    return std::max(0.0, (1 - copyNumber.averageActualBaf) * copyNumber.averageTumorCopyNumber);
}

double majorAlleleCopyNumber(const PurpleCopyNumber &copyNumber){
    return copyNumber.averageTumorCopyNumber - minorAlleleCopyNumber(copyNumber);
}

int deletedGenes(const std::vector<PurpleCopyNumber> &copyNumbers, const std::string &ensemblDataDir){
    struct Exon { int start; int end; };
    std::unordered_map<std::string, std::vector<Exon>> exons;
    std::ifstream exonInput(joinPath(ensemblDataDir, "ensembl_trans_exon_data.csv"));
    if(!exonInput){ throw std::runtime_error("unable to open Ensembl exon data"); }
    std::string line;
    std::getline(exonInput, line);
    while(std::getline(exonInput, line)){
        const auto row = csvFields(line);
        if(row.size() > 10 && row[1] == row[3]){
            exons[row[0]].push_back({std::stoi(row[9]), std::stoi(row[10])});
        }
    }

    std::unordered_map<std::string, std::vector<const PurpleCopyNumber *>> byChromosome;
    for(const auto &copyNumber : copyNumbers){ byChromosome[lp::stripChrPrefix(copyNumber.chromosome)].push_back(&copyNumber); }

    std::ifstream geneInput(joinPath(ensemblDataDir, "ensembl_gene_data.csv"));
    if(!geneInput){ throw std::runtime_error("unable to open Ensembl gene data"); }
    std::getline(geneInput, line);
    int count = 0;
    while(std::getline(geneInput, line)){
        const auto row = csvFields(line);
        if(row.size() < 6){ continue; }
        const std::string chromosome = lp::stripChrPrefix(row[2]);
        if(chromosome == "Y"){ continue; }
        const auto exonIt = exons.find(row[0]);
        const auto copyIt = byChromosome.find(chromosome);
        if(exonIt == exons.end() || copyIt == byChromosome.end()){ continue; }
        double minimum = HUGE_VAL;
        for(const Exon &exon : exonIt->second){
            for(const PurpleCopyNumber *copyNumber : copyIt->second){
                if(copyNumber->start > exon.end){ break; }
                if(copyNumber->end >= exon.start){ minimum = std::min(minimum, copyNumber->averageTumorCopyNumber); }
            }
        }
        if(minimum < 0.5 - EPS){ ++count; }
    }
    return count;
}

double lohPercent(const std::vector<PurpleCopyNumber> &copyNumbers){
    long long totalLohBases = 0;
    long long totalBases = 0;
    for(int chromosome = 1; chromosome <= 22; ++chromosome){
        const std::string name = std::to_string(chromosome);
        std::vector<const PurpleCopyNumber *> values;
        for(const auto &copyNumber : copyNumbers){
            if(lp::stripChrPrefix(copyNumber.chromosome) == name){ values.push_back(&copyNumber); }
        }
        int lohBases = 0;
        int chromosomeBases = 0;
        bool inLoh = false, beforeCentromere = false, hasNonLoh = false;
        int lohStart = 0;
        for(std::size_t i = 0; i < values.size(); ++i){
            const auto &copyNumber = *values[i];
            if(copyNumber.segmentStartSupport == SegmentSupport::TELOMERE){ inLoh = false; beforeCentromere = true; }
            else if(copyNumber.segmentStartSupport == SegmentSupport::CENTROMERE){ beforeCentromere = false; }
            if(beforeCentromere && (chromosome == 13 || chromosome == 14 || chromosome == 15 || chromosome == 21 || chromosome == 22)){ continue; }
            chromosomeBases += copyNumber.end - copyNumber.start + 1;
            const bool isLoh = minorAlleleCopyNumber(copyNumber) < 0.5 - EPS;
            hasNonLoh = hasNonLoh || !isLoh;
            const bool nextIsLoh = i + 1 < values.size() && minorAlleleCopyNumber(*values[i + 1]) < 0.5 - EPS;
            if(!inLoh && isLoh){ inLoh = true; lohStart = copyNumber.start; }
            const bool endOfArm = copyNumber.segmentEndSupport == SegmentSupport::TELOMERE || copyNumber.segmentEndSupport == SegmentSupport::CENTROMERE;
            bool checkEnd = false;
            if(inLoh){
                if(endOfArm){ checkEnd = true; }
                else if(!nextIsLoh){
                    const bool shortInsert = i + 1 < values.size() && values[i + 1]->end - values[i + 1]->start + 1 <= 1000;
                    checkEnd = !shortInsert;
                }
            }
            if(checkEnd){
                const int length = copyNumber.end - lohStart + 1;
                if(length > 1000){ lohBases += length; }
                inLoh = false;
            }
        }
        if(hasNonLoh){ totalLohBases += lohBases; totalBases += chromosomeBases; }
    }
    return totalBases > 0 ? static_cast<double>(totalLohBases) / totalBases : 0;
}

double polyclonalProportion(const std::vector<PurpleCopyNumber> &copyNumbers){
    long long total = 0, polyclonal = 0;
    for(const auto &copyNumber : copyNumbers){
        total += copyNumber.bafCount;
        const double remainder = std::abs(copyNumber.averageTumorCopyNumber - std::round(copyNumber.averageTumorCopyNumber));
        if(remainder > 0.25 + EPS){ polyclonal += copyNumber.bafCount; }
    }
    return total == 0 ? 0 : static_cast<double>(polyclonal) / total;
}

bool wholeGenomeDuplication(const std::vector<PurpleCopyNumber> &copyNumbers){
    int duplicated = 0;
    for(int chromosome = 1; chromosome <= 22; ++chromosome){
        long long count = 0;
        double weighted = 0;
        const std::string name = std::to_string(chromosome);
        for(const auto &copyNumber : copyNumbers){
            if(lp::stripChrPrefix(copyNumber.chromosome) == name){
                count += copyNumber.bafCount;
                weighted += majorAlleleCopyNumber(copyNumber) * copyNumber.bafCount;
            }
        }
        if(count > 0 && weighted / count >= 1.5 - EPS){ ++duplicated; }
    }
    return duplicated >= 11;
}

std::string boolean(bool value){ return value ? "true" : "false"; }

}

SummaryContext buildSummaryContext(const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const std::string &ensemblDataDir){
    SummaryContext result;
    result.deletedGenes = deletedGenes(copyNumbers, ensemblDataDir);
    result.lohPercent = lohPercent(copyNumbers);
    result.polyclonalProportion = polyclonalProportion(copyNumbers);
    result.wholeGenomeDuplication = bestFit.method == "NORMAL" && wholeGenomeDuplication(copyNumbers);
    std::vector<std::string> statuses;
    if(result.deletedGenes > 280){ statuses.push_back("WARN_DELETED_GENES"); }
    if(bestFit.fit.purity < 0.2 - EPS && bestFit.method != "NO_TUMOR"){ statuses.push_back("WARN_LOW_PURITY"); }
    if(inputs.contamination > 0.1 + EPS){ statuses.push_back("FAIL_CONTAMINATION"); }
    if(bestFit.method == "NO_TUMOR"){ statuses.push_back("FAIL_NO_TUMOR"); }
    if(!statuses.empty()){
        result.qcStatus = statuses.front();
        for(std::size_t i = 1; i < statuses.size(); ++i){ result.qcStatus += "," + statuses[i]; }
    }
    return result;
}

void dumpSummaryContext(const InputData &inputs, const BestFit &bestFit,
        const std::vector<PurpleCopyNumber> &copyNumbers, const SummaryContext &summary){
    std::string checkpointStatus = summary.qcStatus;
    for(std::size_t position = 0; (position = checkpointStatus.find(',', position)) != std::string::npos; position += 2){
        checkpointStatus.replace(position, 1, ", ");
    }
    std::vector<lp::CpDump::Row> rows;
    rows.push_back({"", 0,
            std::string(genderName(inputs.amberGender)) + "\tTUMOR\tfalse\t" + bestFit.method + "\t" +
            lp::CpDump::num(bestFit.fit.purity) + "\t" + lp::CpDump::num(bestFit.fit.normFactor) + "\t" +
            lp::CpDump::num(bestFit.fit.ploidy) + "\t" + lp::CpDump::num(bestFit.fit.score) + "\t" +
            lp::CpDump::num(bestFit.fit.diploidProportion) + "\t" + lp::CpDump::num(bestFit.fit.somaticPenalty) +
            "\t[" + checkpointStatus + "]\t" + std::to_string(copyNumbers.size()) + "\t0\t" + std::to_string(summary.deletedGenes) + "\t" +
            lp::CpDump::num(inputs.contamination) + "\t" + std::to_string(inputs.averageTumorDepth) + "\t" +
            genderName(inputs.cobaltGender) + "\t" + genderName(inputs.amberGender) + "\t[NONE]\t" +
            lp::CpDump::num(summary.lohPercent) + "\t0\t" + lp::CpDump::num(summary.polyclonalProportion) + "\t" +
            boolean(summary.wholeGenomeDuplication)});
    lp::CpDump::write("CP-P9-context-qc",
            "gender\trunMode\ttargeted\tmethod\tpurity\tnormFactor\tploidy\tscore\tdiploidProportion\tsomaticPenalty\tqcStatus\tcopyNumberSegments\tunsupportedCopyNumberSegments\tdeletedGenes\tcontamination\tamberMeanDepth\tcobaltGender\tamberGender\tgermlineAberrations\tlohPercent\ttincLevel\tpolyClonalProportion\twholeGenomeDuplication", rows);
}

}
