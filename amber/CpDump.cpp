#include "CpDump.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>

namespace amber {

namespace {
std::string gDir;
}

bool CpDump::enabled(){
    return !gDir.empty();
}

void CpDump::setDir(const std::string &dir){
    gDir = dir;
    if(!gDir.empty()){
        mkdir(gDir.c_str(), 0755);
    }
}

std::string CpDump::num(double value){
    // 17 位有效位數保證 double 可 round-trip；比對端做數值比較，格式不需與 Java 一致
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
}

void CpDump::writeSorted(const std::string &checkpoint, const std::string &header, std::vector<Row> &rows){
    if(!enabled()){
        return;
    }

    // 必須是 stable_sort：Java 端用 List.sort（TimSort，穩定），
    // 同一 (chromosome, position) 的多筆記錄要保留讀入順序
    std::stable_sort(rows.begin(), rows.end(), [](const Row &a, const Row &b){
        if(a.chromosome != b.chromosome){
            return a.chromosome < b.chromosome;
        }
        return a.position < b.position;
    });

    write(checkpoint, header, rows);
}

void CpDump::write(const std::string &checkpoint, const std::string &header, const std::vector<Row> &rows){
    if(!enabled()){
        return;
    }

    const std::string path = gDir + "/" + checkpoint + ".tsv";
    std::ofstream out(path);
    if(!out){
        throw std::runtime_error("cpdump failed to open " + path);
    }

    out << header << '\n';
    for(const Row &row : rows){
        out << row.line << '\n';
    }
}

}
