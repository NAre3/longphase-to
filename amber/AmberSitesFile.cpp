#include "AmberSitesFile.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>

#include <zlib.h>

#include "HumanChromosome.h"

namespace amber {

namespace {

// 讀一行；支援 .gz 與純文字（gzopen 對未壓縮檔亦可讀）
bool readLine(gzFile file, std::string &line){
    line.clear();
    char buffer[65536];

    while(gzgets(file, buffer, sizeof(buffer)) != nullptr){
        line += buffer;
        if(!line.empty() && line.back() == '\n'){
            line.pop_back();
            if(!line.empty() && line.back() == '\r'){
                line.pop_back();
            }
            return true;
        }
    }

    return !line.empty();
}

// 對應 Java 的 line.split("\t", -1)：保留尾端空欄位
std::vector<std::string> splitTab(const std::string &line){
    std::vector<std::string> fields;
    size_t start = 0;

    while(true){
        const size_t tab = line.find('\t', start);
        if(tab == std::string::npos){
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }

    return fields;
}

}

std::vector<AmberSite> loadAmberSites(const std::string &filename){
    gzFile file = gzopen(filename.c_str(), "rb");
    if(file == nullptr){
        throw std::runtime_error("unable to open amber sites file: " + filename);
    }

    std::string header;
    if(!readLine(file, header)){
        gzclose(file);
        throw std::runtime_error("empty amber sites file: " + filename);
    }

    if(header.find("fileformat=VCF") != std::string::npos){
        gzclose(file);
        // Java 端此時改走 loadVcf。本移植階段的輸入是 TSV bundle，遇到 VCF 直接中止，
        // 不做靜默降級——靜默降級會讓兩端跑不同路徑而比對失去意義。
        throw std::runtime_error("VCF amber sites input is not supported by this port: " + filename);
    }

    std::map<std::string, int> fieldIndex;
    const std::vector<std::string> headerFields = splitTab(header);
    for(size_t i = 0; i < headerFields.size(); ++i){
        fieldIndex[headerFields[i]] = static_cast<int>(i);
    }

    auto require = [&](const std::string &name){
        auto iter = fieldIndex.find(name);
        if(iter == fieldIndex.end()){
            gzclose(file);
            throw std::runtime_error("amber sites file missing column: " + name);
        }
        return iter->second;
    };

    const int chrIndex = require("Chromosome");
    const int posIndex = require("Position");
    const int refIndex = require("Ref");
    const int altIndex = require("Alt");
    const int snpCheckIndex = require("SnpCheck");

    // Frequency 為選用欄位；bundle 版本沒有這一欄，Java 端此時以 0 代入
    //
    // **已知限制（2026-09-02 程式碼審查 F-R1）**：此處讀進來的 gnomadFrequency
    // 在下游會被丟棄——PositionEvidence 沒有頻率欄位，NoiseFloor 也把 gnomad 檢查
    // 的 baseline 寫死為 0（NoiseFloor.cpp 的 baselineHetGnomadFrequency）。
    // 因此若換用**帶 Frequency 欄**的 loci 檔，Java 會恢復完整的 gnomad 檢查
    // （PeakGnomadFrequenciesChecker 的兩個 |平均 - baseline| > 0.15 子句），
    // 而本移植不會，且不會報錯——會靜默讓本該被拒絕的 peak 通過。
    //
    // 已知 owner 判定（2026-09-02）：HMF 官方 bundle 的 AmberGermlineSites 不含
    // 此欄且不易變動，故不加防呆拋錯，僅在此註記。改用其他 loci 檔前必須先處理。
    // 完整分析見 research/studies/purple-port-amber-fidelity-v1/CODE-REVIEW-fcb15ca.md
    const bool hasFrequency = fieldIndex.count("Frequency") > 0;
    const int freqIndex = hasFrequency ? fieldIndex["Frequency"] : -1;

    std::vector<AmberSite> sites;
    std::string line;

    while(readLine(file, line)){
        if(line.empty()){
            continue;
        }

        const std::vector<std::string> values = splitTab(line);

        const std::string &chromosome = values[chrIndex];

        // 對應 Java：非 human chromosome 直接跳過（MT、alt contig 等）
        if(!isHumanChromosome(chromosome)){
            continue;
        }

        AmberSite site;
        site.chromosome = chromosome;
        site.position = std::atoi(values[posIndex].c_str());
        site.ref = values[refIndex];
        site.alt = values[altIndex];
        // 對應 Boolean.parseBoolean：只有不分大小寫的 "true" 為真，其餘一律 false
        site.snpCheck = values[snpCheckIndex].size() == 4
                && (values[snpCheckIndex][0] == 't' || values[snpCheckIndex][0] == 'T')
                && (values[snpCheckIndex][1] == 'r' || values[snpCheckIndex][1] == 'R')
                && (values[snpCheckIndex][2] == 'u' || values[snpCheckIndex][2] == 'U')
                && (values[snpCheckIndex][3] == 'e' || values[snpCheckIndex][3] == 'E');
        site.gnomadFrequency = hasFrequency ? std::atof(values[freqIndex].c_str()) : 0.0;

        sites.push_back(site);
    }

    gzclose(file);
    return sites;
}

}
