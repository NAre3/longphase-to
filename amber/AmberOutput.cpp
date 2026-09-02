#include "AmberOutput.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

#include <zlib.h>

namespace amber {

namespace {

// 對應 Java DecimalFormat("0.0000")。
// Java 的 DecimalFormat 預設為 HALF_EVEN；glibc 的 printf %.4f 亦以現行捨入模式
// （round-to-nearest-even）對二進位值做正確捨入。兩者在本研究的資料上是否一致
// 由「與參考輸出逐值比對」實測確認，不以推論代替。
std::string format4(double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);

    // Java 的 DecimalFormat("0.0000") 對 -0.0 會輸出 "-0.0000"；printf 亦然，行為一致。
    return buffer;
}

// AmberBAFFile.toString：normalBAF / normalModifiedBAF 非有限值時寫字面 "0"
// （AmberBAFFile.java:84-86 的 Doubles.isFinite 判斷）。
std::string formatOrZero(double value)
{
    if(!std::isfinite(value)){
        return "0";
    }
    return format4(value);
}

}

double AmberBAF::tumorModifiedBAF() const
{
    return 0.5 + std::fabs(tumorBAF - 0.5);
}

double AmberBAF::normalModifiedBAF() const
{
    return 0.5 + std::fabs(normalBAF - 0.5);
}

void writeAmberBafFile(const std::string &path, const std::vector<AmberBAF> &bafs)
{
    gzFile out = gzopen(path.c_str(), "wb");
    if(out == nullptr){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    // AmberBAFFile.header()
    const std::string header =
            "chromosome\tposition\ttumorBAF\ttumorModifiedBAF\ttumorDepth"
            "\tnormalBAF\tnormalModifiedBAF\tnormalDepth\n";
    if(gzwrite(out, header.data(), static_cast<unsigned>(header.size())) == 0){
        gzclose(out);
        throw std::runtime_error("gzwrite failed: " + path);
    }

    std::string block;
    block.reserve(1 << 20);

    for(const AmberBAF &baf : bafs){
        block += baf.chromosome;
        block += '\t';
        block += std::to_string(baf.position);
        block += '\t';
        block += format4(baf.tumorBAF);
        block += '\t';
        block += format4(baf.tumorModifiedBAF());
        block += '\t';
        block += std::to_string(baf.tumorDepth);
        block += '\t';
        block += formatOrZero(baf.normalBAF);
        block += '\t';
        block += formatOrZero(baf.normalModifiedBAF());
        block += '\t';
        block += std::to_string(baf.normalDepth);
        block += '\n';

        if(block.size() > (1 << 20)){
            if(gzwrite(out, block.data(), static_cast<unsigned>(block.size())) == 0){
                gzclose(out);
                throw std::runtime_error("gzwrite failed: " + path);
            }
            block.clear();
        }
    }

    if(!block.empty()){
        if(gzwrite(out, block.data(), static_cast<unsigned>(block.size())) == 0){
            gzclose(out);
            throw std::runtime_error("gzwrite failed: " + path);
        }
    }

    gzclose(out);
}

void writeAmberQcFile(const std::string &path, double contamination, double consanguinityProportion)
{
    // AmberQC.status()（AmberQC.java:5-16）：
    //   Doubles.greaterThan(contamination, 0.1) → FAIL
    //   Doubles.greaterThan(contamination, 0)   → WARN
    //   否則 PASS
    // Doubles.greaterThan 是 value - reference > 1e-10，不是 >。
    const double epsilon = 1e-10;
    std::string status = "PASS";
    if(contamination - 0.1 > epsilon){
        status = "FAIL";
    }else if(contamination - 0.0 > epsilon){
        status = "WARN";
    }

    std::ofstream out(path);
    if(!out){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    // AmberQCFile.toLines：四行，UniparentalDisomy 在 tumor-only 下為 null → "NONE"
    out << "QCStatus\t" << status << '\n';
    out << "Contamination\t" << format4(contamination) << '\n';
    out << "ConsanguinityProportion\t" << format4(consanguinityProportion) << '\n';
    out << "UniparentalDisomy\tNONE" << '\n';
}


void writeTumorRawFile(const std::string &path, const std::vector<RawTumorRow> &rows)
{
    gzFile out = gzopen(path.c_str(), "wb");
    if(out == nullptr){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    // PositionEvidenceFile.Columns 的順序與名稱（注意是 RefCount/AltCount，不是 refSupport/altSupport）
    const std::string header =
            "Chromosome\tPosition\tRef\tAlt\tReadDepth\tIndelCount"
            "\tRefCount\tAltCount\tBaseQualFiltered\tMapQualFiltered\tSeqTechFiltered\n";
    gzwrite(out, header.data(), static_cast<unsigned>(header.size()));

    std::string block;
    block.reserve(1 << 20);

    for(const RawTumorRow &r : rows){
        block += *r.chromosome;
        block += '\t'; block += std::to_string(r.position);
        block += '\t'; block += r.ref;
        block += '\t'; block += r.alt;
        block += '\t'; block += std::to_string(r.readDepth);
        block += '\t'; block += std::to_string(r.indelCount);
        block += '\t'; block += std::to_string(r.refSupport);
        block += '\t'; block += std::to_string(r.altSupport);
        block += '\t'; block += std::to_string(r.baseQualFiltered);
        block += '\t'; block += std::to_string(r.mapQualFiltered);
        block += '\t'; block += std::to_string(r.seqTechFiltered);
        block += '\n';

        if(block.size() > (1 << 20)){
            gzwrite(out, block.data(), static_cast<unsigned>(block.size()));
            block.clear();
        }
    }

    if(!block.empty()){
        gzwrite(out, block.data(), static_cast<unsigned>(block.size()));
    }

    gzclose(out);
}

void writeVersionFile(const std::string &path)
{
    std::ofstream out(path);
    if(!out){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    // Java 端的 amber.version 為 "version=4.3" + "build.date=<jar 建置時間>"。
    // 此處刻意不冒用該版本字串，且 build.date 本就無法重現（那是 jar 的建置時間）。
    out << "version=amber-4.3-cpp-port\n";
    out << "reproduces=hmftools amber-v4.3\n";
    out << "implementation=LongPhase-TO C++ port\n";
}

}
