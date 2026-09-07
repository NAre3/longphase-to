#ifndef COBALT_GCBUCKETS_H
#define COBALT_GCBUCKETS_H

#include <cmath>
#include <vector>

namespace cobalt {

// GCPail.java
class GcPail
{
public:
    explicit GcPail(int gc) : mGc(gc) {}

    // bucketIndex(gc) = (int) Math.round(gc * 100)
    // Math.round 為 half-up 且回傳 long；gc ∈ [0,1] 為非負，與 std::lround 等價。
    static int bucketIndex(double gcContent) { return static_cast<int>(std::lround(gcContent * 100)); }

    int gc() const { return mGc; }
    long readingCount() const { return static_cast<long>(mValues.size()); }
    void addReading(double value) { mValues.push_back(value); }

    // median() = getN() == 0 ? 0 : getPercentile(50)   ← 特例回 0，不是 NaN
    double median() const;

private:
    int mGc;
    std::vector<double> mValues;
};

class GcPailsList
{
public:
    GcPailsList();
    GcPail &getGcPail(double gcContent) { return mBuckets[static_cast<std::size_t>(GcPail::bucketIndex(gcContent))]; }
    const std::vector<GcPail> &buckets() const { return mBuckets; }
    std::vector<GcPail> &buckets() { return mBuckets; }

private:
    std::vector<GcPail> mBuckets;
};

// GcBucketStatistics.java。**保留原始的 off-by-one**（見 behaviour-contract.md §2.1）。
class GcBucketStatistics
{
public:
    GcBucketStatistics(const GcPailsList &pails, int minAllowedGc, int maxAllowedGc);

    // 直接取 MeanDepths[i]，不經 isAllowed（CP-C8 的 dump 用）
    double medianReadDepthForBucket(int i) const { return mMeanDepths[static_cast<std::size_t>(i)]; }

    // medianReadDepth(pail)：isAllowed 檢查的是 pail 自己，讀的是 MeanDepths[pail.mGC]
    double medianReadDepth(const GcPail *pail) const;

    bool isAllowed(const GcPail *pail) const
    {
        return pail != nullptr && pail->gc() > mMinAllowedGc && pail->gc() <= mMaxAllowedGc;
    }

private:
    int mMinAllowedGc;
    int mMaxAllowedGc;
    std::vector<double> mMeanDepths;
};

}

#endif
