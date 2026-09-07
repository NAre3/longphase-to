#include "CobaltOutput.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <zlib.h>

#include "../common/Segmentation.h"

namespace cobalt {

namespace {

// java.text.DecimalFormat("#.####") 的輸出：
//   - 最多 4 位小數，**去掉尾隨零**（-1.0 -> "-1"、0.0 -> "0"、0.70040 -> "0.7004"）
//   - 預設 RoundingMode.HALF_EVEN（**與 Math.round 的 half-up 不同**）
std::string decimalFormat4(double v)
{
    if(std::isnan(v)){ return "NaN"; }

    // 以 4 位小數的 half-even 捨入；printf 的 %f 在 IEEE-754 下即為 half-even
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4f", v);
    std::string s(buf);

    // 去掉尾隨零與可能剩下的小數點
    if(s.find('.') != std::string::npos)
    {
        while(!s.empty() && s.back() == '0'){ s.pop_back(); }
        if(!s.empty() && s.back() == '.'){ s.pop_back(); }
    }
    // DecimalFormat 不會輸出 "-0"
    if(s == "-0"){ s = "0"; }
    return s;
}

}

std::vector<CobaltRatio> collateResults(const std::vector<BamRatio> &tumorRatios)
{
    // ResultsCollator.tumorOnlyCobaltRatio：四個 reference 欄位固定 -1.0（G14）
    std::vector<CobaltRatio> out;
    out.reserve(tumorRatios.size());
    for(const BamRatio &r : tumorRatios)
    {
        CobaltRatio c;
        c.chromosome = "chr" + r.chromosomeShort;               // versionedChromosome(V38)
        c.position = r.position;
        c.referenceReadDepth = -1.0;
        c.referenceGCRatio = -1.0;
        c.referenceGCContent = -1.0;
        c.referenceGCDiploidRatio = -1.0;
        c.tumorReadDepth = r.readDepth();
        c.tumorGCRatio = r.ratio();
        c.tumorGcContent = r.gcContent();
        out.push_back(std::move(c));
    }
    return out;
}

void writeCobaltRatioFile(const std::string &path, const std::vector<CobaltRatio> &ratios)
{
    gzFile fp = gzopen(path.c_str(), "wb");
    if(fp == nullptr){ throw std::runtime_error("cannot open for write: " + path); }

    // 檔案欄序由 lambda$write$1 的 Row.set 呼叫順序讀出（位元組碼），**與物件欄序不同**
    const char *header = "chromosome\tposition\treferenceReadDepth\ttumorReadDepth\t"
                         "referenceGCRatio\ttumorGCRatio\treferenceGCDiploidRatio\t"
                         "referenceGCContent\ttumorGCContent\n";
    gzwrite(fp, header, static_cast<unsigned>(std::string(header).size()));

    std::string buf;
    buf.reserve(1 << 20);
    for(const CobaltRatio &r : ratios)
    {
        buf += r.chromosome; buf += '\t';
        buf += std::to_string(r.position); buf += '\t';
        buf += decimalFormat4(r.referenceReadDepth); buf += '\t';
        buf += decimalFormat4(r.tumorReadDepth); buf += '\t';
        buf += decimalFormat4(r.referenceGCRatio); buf += '\t';
        buf += decimalFormat4(r.tumorGCRatio); buf += '\t';
        buf += decimalFormat4(r.referenceGCDiploidRatio); buf += '\t';
        buf += decimalFormat4(r.referenceGCContent); buf += '\t';
        buf += decimalFormat4(r.tumorGcContent); buf += '\n';
        if(buf.size() > (1u << 19))
        {
            gzwrite(fp, buf.data(), static_cast<unsigned>(buf.size()));
            buf.clear();
        }
    }
    if(!buf.empty()){ gzwrite(fp, buf.data(), static_cast<unsigned>(buf.size())); }
    gzclose(fp);
}

void writeGcMedianFile(const std::string &path, double sampleMean, double sampleMedian,
                       const GcBucketStatistics &stats)
{
    FILE *f = std::fopen(path.c_str(), "w");
    if(f == nullptr){ throw std::runtime_error("cannot open for write: " + path); }
    // String.format("%.2f") 為 half-up（Java Formatter 的 RoundingMode.HALF_UP）
    std::fprintf(f, "#sampleMean\tsampleMedian\n");
    std::fprintf(f, "%.2f\t%.2f\n", sampleMean, sampleMedian);
    std::fprintf(f, "#gcBucket\tmedian\n");
    for(int i = 0; i < 101; ++i)
    {
        const double m = stats.medianReadDepthForBucket(i);
        if(m > 0){ std::fprintf(f, "%d\t%.2f\n", i, m); }        // bucketToMedianReadDepth: MeanDepths[i] > 0
    }
    std::fclose(f);
}

}
