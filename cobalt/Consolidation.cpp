#include "Consolidation.h"

#include <cmath>
#include <stdexcept>

#include "CobaltConstants.h"
#include "../common/Segmentation.h"

namespace cobalt {

namespace {

// roundDownToWindowBoundary(p) = (int)(Math.floor(p / 1000) * 1000) + 1
int roundDownToWindowBoundary(double p)
{
    return static_cast<int>(std::floor(p / WINDOW_SIZE) * WINDOW_SIZE) + 1;
}

// commons-math3 Mean.evaluate：兩趟修正式平均（與 BamRatio.cpp 的同一形式）。
// Sum.evaluate 為單純循序累加，無補償。
double meanEvaluate(const std::vector<double> &v)
{
    if(v.empty()){ return std::nan(""); }
    double sum = 0.0;
    for(double x : v){ sum += x; }
    const double sampleSize = static_cast<double>(v.size());
    const double xbar = sum / sampleSize;
    double correction = 0.0;
    for(double x : v){ correction += (x - xbar); }
    return xbar + correction / sampleSize;
}

}

int calcConsolidationCount(double medianReadDepth)
{
    const double c = 80.0 / medianReadDepth;
    if(c < 10.0){ return 1; }                                   // 嚴格小於，即 median > 8
    if(c >= 1000.0){ return 1000; }
    const double roundBy = std::pow(10.0, std::floor(std::log10(c)));
    // Math.round 為 half-up 回傳 long；(int) 截斷發生在乘 roundBy 之後
    return static_cast<int>(std::llround(c / roundBy) * roundBy);
}

ConsolidatorChoice chooseConsolidator(double medianReadDepth)
{
    if(std::isnan(medianReadDepth)){ return ConsolidatorChoice{"NoOpConsolidator", 1}; }
    const int count = calcConsolidationCount(medianReadDepth);
    if(count == 1){ return ConsolidatorChoice{"NoOpConsolidator", 1}; }
    return ConsolidatorChoice{"LowCoverageConsolidator", count};
}

std::vector<LowCovBucket> consolidateIntoBuckets(const std::vector<int> &windowPositions, int consolidationCount)
{
    std::vector<LowCovBucket> buckets;
    if(windowPositions.empty()){ return buckets; }

    // Validate.isTrue(Comparators.isInStrictOrder(...))
    for(std::size_t i = 1; i < windowPositions.size(); ++i)
    {
        if(windowPositions[i] <= windowPositions[i - 1])
        {
            throw std::runtime_error("window positions not in strict order");
        }
    }

    int windowCount = 0;
    int bucketStart = windowPositions[0];

    for(std::size_t i = 0; i < windowPositions.size(); ++i)
    {
        const int position = windowPositions[i];

        // 兩個 if 依序判斷，不是 else if
        if((position - bucketStart) >= MAX_SPARSE_CONSOLIDATE_DISTANCE)
        {
            if(i > 0)
            {
                const int lastPosition = windowPositions[i - 1];
                const int bucketEnd = lastPosition + WINDOW_SIZE;
                const int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
                buckets.push_back(LowCovBucket{bucketStart, bucketEnd, bucketPos});
            }
            bucketStart = position;
            windowCount = 0;
        }

        if(windowCount == consolidationCount)
        {
            const int lastPosition = windowPositions[i - 1];
            const int bucketEnd = roundDownToWindowBoundary((lastPosition + position) * 0.5);
            const int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
            buckets.push_back(LowCovBucket{bucketStart, bucketEnd, bucketPos});
            bucketStart = bucketEnd + WINDOW_SIZE;
            windowCount = 0;
        }

        ++windowCount;                                          // 在迴圈尾端
    }

    if(windowCount > 0)
    {
        const int bucketEnd = windowPositions.back() + WINDOW_SIZE;
        const int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
        buckets.push_back(LowCovBucket{bucketStart, bucketEnd, bucketPos});
    }

    return buckets;
}

void applyLowCoverageConsolidation(std::vector<BamRatio> &ratios, int consolidationCount)
{
    // ---- 依染色體切段（ratios 已依 ordinal、position 排序）----
    std::vector<std::pair<std::size_t, std::size_t>> spans;     // [begin, end)
    for(std::size_t i = 0; i < ratios.size();)
    {
        std::size_t j = i;
        while(j < ratios.size() && ratios[j].chromosomeShort == ratios[i].chromosomeShort){ ++j; }
        spans.emplace_back(i, j);
        i = j;
    }

    std::size_t consolidatedTotal = 0;
    std::vector<std::vector<BamRatio>> perChromosome(spans.size());

    for(std::size_t s = 0; s < spans.size(); ++s)
    {
        const std::size_t b = spans[s].first, e = spans[s].second;

        // createBoundariesIfNecessary：非遮蔽（ratio() >= 0，非嚴格）的 position
        std::vector<int> nonMasked;
        for(std::size_t i = b; i < e; ++i)
        {
            if(ratios[i].ratio() >= 0){ nonMasked.push_back(ratios[i].position); }
        }
        const std::vector<LowCovBucket> boundaries = consolidateIntoBuckets(nonMasked, consolidationCount);
        if(boundaries.empty()){ continue; }                     // Java：bucketItr 無元素 -> log error + continue

        // populateLowCoverageRatio
        std::vector<BamRatio> &out = perChromosome[s];
        std::size_t bi = 0;
        bool haveBucket = true;
        std::vector<double> bucketRatios, bucketGcs;

        auto emit = [&](const LowCovBucket &bk)
        {
            // new BamRatio(chromosome, bucket.BucketPosition, ratios.getMean(), gcs.getMean())
            // 走第二個建構子：ratio 初值 = readDepth = bucketRatios.getMean()
            BamRatio r;
            r.chromosomeShort = ratios[b].chromosomeShort;
            r.position = bk.bucketPosition;
            r.setFromConsolidated(meanEvaluate(bucketRatios), meanEvaluate(bucketGcs));
            out.push_back(std::move(r));
        };

        for(std::size_t i = b; i < e; ++i)
        {
            if(!haveBucket){ continue; }
            if(!(ratios[i].ratio() >= 0)){ continue; }
            if(ratios[i].position > boundaries[bi].endPosition)
            {
                emit(boundaries[bi]);
                if(bi + 1 < boundaries.size())
                {
                    ++bi;
                    bucketRatios.clear(); bucketGcs.clear();
                    // 觸發前進的那一列會被加進新 bucket
                    bucketRatios.push_back(ratios[i].ratio());
                    bucketGcs.push_back(ratios[i].gcContent());
                }
                else
                {
                    haveBucket = false;
                }
            }
            else
            {
                bucketRatios.push_back(ratios[i].ratio());
                bucketGcs.push_back(ratios[i].gcContent());
            }
        }
        if(haveBucket){ emit(boundaries[bi]); }
        consolidatedTotal += out.size();
    }

    // ---- BamRatios.consolidate 的回填 ----
    if(consolidatedTotal == ratios.size()){ return; }            // size 相等 -> 提前返回

    for(std::size_t s = 0; s < spans.size(); ++s)
    {
        const std::size_t b = spans[s].first, e = spans[s].second;
        std::size_t oi = b;
        for(const BamRatio &c : perChromosome[s])
        {
            bool seeking = true;
            while(seeking && oi < e)
            {
                BamRatio &o = ratios[oi++];
                if(o.position == c.position){ o.overrideRatio(c.ratio()); seeking = false; }
                else                        { o.overrideRatio(-1.0); }
            }
        }
        while(oi < e){ ratios[oi++].overrideRatio(-1.0); }        // 尾端全 mask
    }
}

}
