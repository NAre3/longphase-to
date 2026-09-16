#include "PurpleInput.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <zlib.h>

#include "../common/CpDump.h"
#include "../common/HumanChromosome.h"

namespace purple {
namespace {

std::vector<std::string> split(const std::string &line){
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while(true){
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == std::string::npos ? end : end - begin));
        if(end == std::string::npos){ break; }
        begin = end + 1;
    }
    return fields;
}

std::unordered_map<std::string, std::size_t> columns(const std::string &header){
    std::unordered_map<std::string, std::size_t> result;
    const auto values = split(header);
    for(std::size_t i = 0; i < values.size(); ++i){ result.emplace(values[i], i); }
    return result;
}

template<typename Consumer>
void forEachGzipLine(const std::string &path, Consumer consumer){
    gzFile input = gzopen(path.c_str(), "rb");
    if(input == nullptr){ throw std::runtime_error("unable to open gzip input: " + path); }
    std::string line;
    char buffer[65536];
    while(gzgets(input, buffer, sizeof(buffer)) != nullptr){
        line.append(buffer);
        if(!line.empty() && line.back() == '\n'){
            line.pop_back();
            if(!line.empty() && line.back() == '\r'){ line.pop_back(); }
            consumer(line);
            line.clear();
        }
    }
    if(!line.empty()){ consumer(line); }
    const int status = gzclose(input);
    if(status != Z_OK){ throw std::runtime_error("error reading gzip input: " + path); }
}

std::size_t gzipDataLineCount(const std::string &path){
    std::size_t lines = 0;
    forEachGzipLine(path, [&](const std::string &){ ++lines; });
    return lines == 0 ? 0 : lines - 1;
}

bool isChr(const std::string &value, const char *name){
    return lp::stripChrPrefix(value) == name;
}

double genderAdjusted(Gender gender, const std::string &chromosome, double value){
    if(isChr(chromosome, "X")){ return gender == Gender::FEMALE ? 1.0 : 0.5; }
    if(isChr(chromosome, "Y")){ return gender == Gender::FEMALE ? 0.0 : 0.5; }
    return value;
}

std::string joinPath(const std::string &directory, const std::string &name){
    return directory + (directory.empty() || directory.back() == '/' ? "" : "/") + name;
}

std::string boolean(bool value){ return value ? "true" : "false"; }

}

std::vector<AmberBaf> readAmberBafs(const std::string &path){
    std::vector<AmberBaf> result;
    result.reserve(gzipDataLineCount(path));
    std::unordered_map<std::string, std::size_t> col;
    bool header = true;
    forEachGzipLine(path, [&](const std::string &line){
        if(header){ col = columns(line); header = false; return; }
        const auto row = split(line);
        AmberBaf baf;
        baf.chromosome = row.at(col.at("chromosome"));
        baf.position = std::stoi(row.at(col.at("position")));
        baf.tumorBaf = std::stod(row.at(col.at("tumorBAF")));
        baf.tumorDepth = std::stoi(row.at(col.at("tumorDepth")));
        baf.normalBaf = std::stod(row.at(col.at("normalBAF")));
        baf.normalDepth = std::stoi(row.at(col.at("normalDepth")));
        result.push_back(std::move(baf));
    });
    if(header){ throw std::runtime_error("empty Amber BAF input: " + path); }
    return result;
}

std::vector<CobaltRatio> readCobaltRatios(const std::string &path, Gender gender){
    std::vector<CobaltRatio> result;
    result.reserve(gzipDataLineCount(path));
    std::unordered_map<std::string, std::size_t> col;
    bool header = true;
    forEachGzipLine(path, [&](const std::string &line){
        if(header){ col = columns(line); header = false; return; }
        const auto row = split(line);
        CobaltRatio ratio;
        ratio.chromosome = row.at(col.at("chromosome"));
        ratio.position = std::stoi(row.at(col.at("position")));
        ratio.referenceReadDepth = std::stod(row.at(col.at("referenceReadDepth")));
        double referenceGcRatio = std::stod(row.at(col.at("referenceGCRatio")));
        double referenceGcDiploidRatio = std::stod(row.at(col.at("referenceGCDiploidRatio")));
        if(ratio.referenceReadDepth == -1){
            referenceGcRatio = 1;
            referenceGcDiploidRatio = 1;
        }
        ratio.referenceGcRatio = genderAdjusted(gender, ratio.chromosome, referenceGcRatio);
        ratio.referenceGcDiploidRatio = genderAdjusted(gender, ratio.chromosome, referenceGcDiploidRatio);
        ratio.referenceGcContent = std::stod(row.at(col.at("referenceGCContent")));
        ratio.tumorReadDepth = std::stod(row.at(col.at("tumorReadDepth")));
        ratio.tumorGcRatio = std::stod(row.at(col.at("tumorGCRatio")));
        ratio.tumorGcContent = std::stod(row.at(col.at("tumorGCContent")));
        result.push_back(std::move(ratio));
    });
    if(header){ throw std::runtime_error("empty Cobalt ratio input: " + path); }
    return result;
}

std::vector<PcfPosition> readPcfPositions(const std::string &path, PcfSource source){
    std::ifstream input(path);
    if(!input){ throw std::runtime_error("unable to open PCF input: " + path); }
    std::string line;
    if(!std::getline(input, line)){ throw std::runtime_error("empty PCF input: " + path); }
    const bool oldFormat = line.rfind("sampleID", 0) == 0;
    std::unordered_map<std::string, std::vector<PcfPosition>> byChromosome;
    std::vector<std::string> chromosomeOrder;
    std::string currentChromosome;
    int minPosition = 1;
    while(std::getline(input, line)){
        const auto row = split(line);
        const std::string chromosome = row.at(oldFormat ? 1 : 0);
        if(!lp::isHumanChromosome(chromosome)){ continue; }
        if(chromosome != currentChromosome){
            currentChromosome = chromosome;
            minPosition = 1;
            chromosomeOrder.push_back(chromosome);
        }
        const int rawStart = std::stoi(row.at(oldFormat ? 3 : 1));
        const int rawEnd = std::stoi(row.at(oldFormat ? 4 : 2));
        const int start = oldFormat ? ((rawStart - 1) / 1000) * 1000 + 1 : rawStart;
        const int end = oldFormat ? ((rawEnd - 1) / 1000) * 1000 + 1001 : rawEnd + 1;
        auto &positions = byChromosome[chromosome];
        if(!positions.empty()){ positions.back().maxPosition = start; }
        positions.push_back(PcfPosition{source, chromosome, start, minPosition, start});
        minPosition = end;
        positions.push_back(PcfPosition{source, chromosome, end, end, end});
    }

    std::vector<PcfPosition> result;
    for(const std::string &chromosome : chromosomeOrder){
        auto &positions = byChromosome.at(chromosome);
        std::stable_sort(positions.begin(), positions.end(), [](const PcfPosition &a, const PcfPosition &b){
            return a.position < b.position;
        });
        for(std::size_t i = 0; i + 1 < positions.size(); ){
            auto &current = positions[i];
            auto &next = positions[i + 1];
            if(current.position == next.position){
                current.minPosition = std::max(current.minPosition, next.minPosition);
                current.maxPosition = std::min(current.maxPosition, next.maxPosition);
                positions.erase(positions.begin() + static_cast<std::ptrdiff_t>(i + 1));
            }else{
                current.maxPosition = std::min(current.maxPosition, next.position);
                next.minPosition = std::max(next.minPosition, current.position);
                ++i;
            }
        }
        result.insert(result.end(), positions.begin(), positions.end());
    }
    return result;
}

double readAmberContamination(const std::string &path){
    std::ifstream input(path);
    if(!input){ throw std::runtime_error("unable to open Amber QC input: " + path); }
    std::string line;
    while(std::getline(input, line)){
        const auto row = split(line);
        if(row.size() >= 2 && row[0] == "Contamination"){ return std::stod(row[1]); }
    }
    throw std::runtime_error("Amber QC has no Contamination field: " + path);
}

Gender determineAmberGender38(const std::vector<AmberBaf> &bafs){
    std::size_t xPoints = 0;
    for(const AmberBaf &baf : bafs){
        if(isChr(baf.chromosome, "X") && baf.position >= 2781479 && baf.position <= 156030895){ ++xPoints; }
    }
    if(xPoints < 3){ return Gender::MALE; }
    return static_cast<double>(xPoints) / static_cast<double>(bafs.size()) > 0.01 ? Gender::FEMALE : Gender::MALE;
}

int averageTumorDepth(const std::vector<AmberBaf> &bafs){
    std::int64_t total = 0;
    std::size_t count = 0;
    for(const AmberBaf &baf : bafs){ if(baf.tumorDepth > 0){ total += baf.tumorDepth; ++count; } }
    return count == 0 ? 100 : static_cast<int>(std::floor(static_cast<double>(total) / count + 0.5));
}

InputData loadTumorOnlyInputs(const std::string &sampleId, const std::string &amberDir, const std::string &cobaltDir){
    InputData data;
    data.sampleId = sampleId;
    data.bafs = readAmberBafs(joinPath(amberDir, sampleId + ".amber.baf.tsv.gz"));
    data.contamination = readAmberContamination(joinPath(amberDir, sampleId + ".amber.qc"));
    data.averageTumorDepth = averageTumorDepth(data.bafs);
    data.amberGender = determineAmberGender38(data.bafs);
    data.cobaltGender = data.amberGender;
    data.amberPcf = readPcfPositions(joinPath(amberDir, sampleId + ".amber.baf.pcf"), PcfSource::TUMOR_BAF);
    data.cobaltTumorPcf = readPcfPositions(joinPath(cobaltDir, sampleId + ".cobalt.ratio.pcf"), PcfSource::TUMOR_RATIO);
    data.ratios = readCobaltRatios(joinPath(cobaltDir, sampleId + ".cobalt.ratio.tsv.gz"), data.amberGender);
    return data;
}

void dumpInputCheckpoint(const InputData &data){
    using lp::CpDump;
    std::vector<CpDump::Row> rows;
    const std::string contamination = data.contamination == 0.0 ? "0.0" : CpDump::num(data.contamination);
    rows.push_back({"", 0, data.sampleId + "\tnull\t" + genderName(data.amberGender) + "\t" + genderName(data.cobaltGender) + "\t" +
            std::to_string(data.averageTumorDepth) + "\t" + contamination + "\t" + std::to_string(data.bafs.size()) + "\t" +
            std::to_string(data.ratios.size()) + "\t" + std::to_string(data.amberPcf.size()) + "\t" +
            std::to_string(data.cobaltTumorPcf.size()) + "\t0"});
    CpDump::write("CP-P1-summary", "sample\treference\tamberGender\tcobaltGender\taverageTumorDepth\tcontamination\tbafCount\tratioCount\tamberPcfCount\tcobaltTumorPcfCount\tcobaltReferencePcfCount", rows);

    rows.clear(); rows.reserve(data.bafs.size());
    for(const auto &x : data.bafs){ rows.push_back({x.chromosome, x.position, x.chromosome + "\t" + std::to_string(x.position) + "\t" + CpDump::num(x.tumorBaf) + "\t" + std::to_string(x.tumorDepth) + "\t" + CpDump::num(x.normalBaf) + "\t" + std::to_string(x.normalDepth)}); }
    CpDump::write("CP-P1-amber-baf", "chromosome\tposition\ttumorBAF\ttumorDepth\tnormalBAF\tnormalDepth", rows);

    rows.clear(); rows.reserve(data.ratios.size());
    for(const auto &x : data.ratios){ rows.push_back({x.chromosome, x.position, x.chromosome + "\t" + std::to_string(x.position) + "\t" + CpDump::num(x.referenceReadDepth) + "\t" + CpDump::num(x.referenceGcRatio) + "\t" + CpDump::num(x.referenceGcContent) + "\t" + CpDump::num(x.referenceGcDiploidRatio) + "\t" + CpDump::num(x.tumorReadDepth) + "\t" + CpDump::num(x.tumorGcRatio) + "\t" + CpDump::num(x.tumorGcContent)}); }
    CpDump::write("CP-P1-cobalt-ratio", "chromosome\tposition\treferenceReadDepth\treferenceGCRatio\treferenceGcContent\treferenceGCDiploidRatio\ttumorReadDepth\ttumorGCRatio\ttumorGcContent", rows);

    const auto dumpPcf = [&](const char *name, const std::vector<PcfPosition> &positions){
        std::vector<CpDump::Row> pcfRows; pcfRows.reserve(positions.size());
        for(const auto &x : positions){ pcfRows.push_back({x.chromosome, x.position, std::string(pcfSourceName(x.source)) + "\t" + x.chromosome + "\t" + std::to_string(x.position) + "\t" + std::to_string(x.minPosition) + "\t" + std::to_string(x.maxPosition) + "\t" + boolean(x.isSegmentStart()) + "\t" + boolean(x.isSegmentEnd())}); }
        CpDump::write(name, "source\tchromosome\tposition\tminPosition\tmaxPosition\tisStart\tisEnd", pcfRows);
    };
    dumpPcf("CP-P1-amber-pcf", data.amberPcf);
    dumpPcf("CP-P1-cobalt-tumor-pcf", data.cobaltTumorPcf);
    dumpPcf("CP-P1-cobalt-reference-pcf", {});
}

}
