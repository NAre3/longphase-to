#ifndef COBALT_BAMRATIO_H
#define COBALT_BAMRATIO_H

#include <cmath>
#include <string>
#include <vector>

#include "CobaltWindow.h"

namespace cobalt {

// BamRatio.java
class BamRatio
{
public:
    std::string chromosomeShort;
    int position = 0;

    // CobaltWindow.toBamRatio -> BamRatio(chromosome, position, readDepth, gcContent, included)
    static BamRatio fromWindow(const CobaltWindow &w);

    void normaliseForGc(double medianReadDepthForGcBucket) { normalise(medianReadDepthForGcBucket); }
    void applyEnrichment(double enrichment)                { normalise(enrichment); }
    void normaliseByMean(double mean)                      { normalise(mean); }

    double readDepth() const { return mReadDepth; }
    double gcContent() const { return mGcContent; }
    double ratio() const     { return mIncluded ? mRatio : -1.0; }   // getter 已含遮蔽
    bool included() const    { return mIncluded; }
    double diploidAdjustedRatio() const { return mDiploidAdjustedRatio; }

    void overrideRatio(double r) { mRatio = r; if(mRatio > 0){ mIncluded = true; } }

    // BamRatio(chromosome, position, readDepth, gcContent) -> this(…, readDepth, readDepth, gcContent)
    // 第二個建構子：Ratio 初值 = readDepth，Included = true。LowCoverageConsolidator 用。
    void setFromConsolidated(double readDepth, double gcContent)
    {
        mReadDepth = readDepth; mRatio = readDepth; mGcContent = gcContent; mIncluded = true;
    }

private:
    void normalise(double factor);

    double mReadDepth = 0.0;
    double mRatio = 0.0;
    double mDiploidAdjustedRatio = -1.0;
    double mGcContent = 0.0;
    bool mIncluded = true;
};

// ReadDepthStatisticsNormaliser.java
class ReadDepthStatisticsNormaliser
{
public:
    void recordValue(const BamRatio &r);
    void dataCollectionFinished();
    void normalise(BamRatio &r) const;

    double readDepthMean() const   { return mReadDepthMean; }
    double readDepthMedian() const { return mReadDepthMedian; }
    long sampleCount() const       { return static_cast<long>(mValues.size()); }

    // §5.2 見證欄位 #1：recordValue 是否會把這個 window 加進統計。
    // 述詞逐字對應原碼的三個 early return，不可改寫為 gc >= MIN && gc <= MAX（NaN 上相反）。
    static bool includedInStats(const BamRatio &r);

private:
    std::vector<double> mValues;
    double mRatio = -1.0;
    double mReadDepthMean = -1.0;
    double mReadDepthMedian = -1.0;
};

}

#endif
