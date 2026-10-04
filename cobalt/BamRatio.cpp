#include "BamRatio.h"

#include <algorithm>

#include "CobaltConstants.h"
#include "Percentile.h"

namespace cobalt {

namespace {

// Doubles.equal(a, b) = Math.abs(a - b) < 1.0E-10（嚴格小於）
// Doubles.isZero(x)   = Doubles.equal(x, 0.0)
// **注意這不是 x == 0.0**，帶 1e-10 的 epsilon。
inline bool doublesIsZero(double x){ return std::fabs(x - 0.0) < DOUBLES_EPSILON; }

bool isAutosome(const std::string &shortName)
{
    if(shortName == "X" || shortName == "Y"){ return false; }
    for(char c : shortName){ if(c < '0' || c > '9'){ return false; } }
    return !shortName.empty();
}

// commons-math3 Sum.evaluate：單純循序累加，**無補償**
// （與 EXP-C005 的 DoubleStream.average 的 Kahan 補償不同，不可混用）
double sumEvaluate(const std::vector<double> &v)
{
    double sum = 0.0;
    for(double x : v){ sum += x; }
    return sum;
}

// commons-math3 Mean.evaluate(double[], int, int)：兩趟修正式平均
double meanEvaluate(const std::vector<double> &v)
{
    if(v.empty()){ return std::nan(""); }
    const double sampleSize = static_cast<double>(v.size());
    const double xbar = sumEvaluate(v) / sampleSize;
    double correction = 0.0;
    for(double x : v){ correction += (x - xbar); }
    return xbar + correction / sampleSize;
}

}

BamRatio BamRatio::fromWindow(const CobaltWindow &w)
{
    // BamRatio(chromosome, position, readDepth, gcContent, included) —— 第四個建構子
    BamRatio r;
    r.chromosomeShort = w.chromosomeShort;
    r.position = w.position;
    r.mReadDepth = std::isfinite(w.readDepth) ? w.readDepth : -1.0;
    r.mRatio = r.mReadDepth;                                    // Ratio 初值 = mReadDepth
    r.mGcContent = std::isfinite(w.gcContent) ? w.gcContent : -1.0;
    r.mIncluded = !w.isInExcludedRegion && w.isInTargetRegion;  // CobaltWindow.include()
    if(!r.mIncluded){ r.mRatio = -1.0; }
    return r;
}

void BamRatio::normalise(double factor)
{
    // 四個分支的順序不可調換
    if(std::isnan(factor)){ mIncluded = false; mRatio = -1.0; return; }
    if(doublesIsZero(mRatio)){ return; }                        // 帶 1e-10 epsilon
    if(factor <= 0 || !mIncluded || mRatio < 0)
    {
        mIncluded = false;
        mRatio = -1.0;
    }
    else
    {
        mRatio /= factor;
    }
}

bool ReadDepthStatisticsNormaliser::includedInStats(const BamRatio &r)
{
    // 逐字對應三個 early return；否定形式在 NaN 上與 >=/<= 相反，不可改寫
    if(!isAutosome(r.chromosomeShort)){ return false; }
    if(r.gcContent() < DEFAULT_GC_RATIO_MIN){ return false; }
    if(r.gcContent() > DEFAULT_GC_RATIO_MAX){ return false; }
    return r.ratio() > 0.0;
}

void ReadDepthStatisticsNormaliser::recordValue(const BamRatio &r)
{
    if(!isAutosome(r.chromosomeShort)){ return; }
    if(r.gcContent() < DEFAULT_GC_RATIO_MIN){ return; }
    if(r.gcContent() > DEFAULT_GC_RATIO_MAX){ return; }
    if(r.ratio() > 0.0){ mValues.push_back(r.readDepth()); }    // 加的是 readDepth，不是 ratio
}

void ReadDepthStatisticsNormaliser::dataCollectionFinished()
{
    mReadDepthMean = meanEvaluate(mValues);
    std::vector<double> copy = mValues;                          // percentile 會就地重排
    mReadDepthMedian = mValues.empty() ? std::nan("") : percentile(copy, 50.0);
    mRatio = mReadDepthMean / mReadDepthMedian;
}

void ReadDepthStatisticsNormaliser::normalise(BamRatio &r) const
{
    r.normaliseByMean(mRatio);
}

}
