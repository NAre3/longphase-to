#include "GcBuckets.h"

#include "Percentile.h"

namespace {

// java.util.DoubleStream.average() -> DoubleSummaryStatistics
// accept() 逐項呼叫 sumWithCompensation（Kahan 補償），getSum() 回傳 sum - sumCompensation，
// getAverage() 再除以 count。**屬 D1 的「加總的形式」＝演算法邏輯，必須逐字重現**——
// 與 commons-math Mean 的兩趟修正式平均同理。
// 實測：單純 (a+b+c)/3 只在 44 格中的 42 格重現 Java；本形式 44/44 全中。
double doubleStreamAverage(const double *values, int count)
{
    double sum = 0.0;
    double sumCompensation = 0.0;
    for(int i = 0; i < count; ++i)
    {
        const double tmp = values[i] - sumCompensation;
        const double velvel = sum + tmp;
        sumCompensation = (velvel - sum) - tmp;
        sum = velvel;
    }
    // getSum 的 NaN/Infinite 退回 simpleSum 分支在此不可能觸發（三個有限值）
    return (sum - sumCompensation) / count;
}

}

namespace cobalt {

double GcPail::median() const
{
    if(mValues.empty()){ return 0; }          // getN() == 0 → 0（不是 NaN）
    std::vector<double> copy = mValues;       // percentile 會就地重排
    return percentile(copy, 50.0);
}

GcPailsList::GcPailsList()
{
    mBuckets.reserve(101);
    for(int i = 0; i < 101; ++i){ mBuckets.emplace_back(i); }
}

GcBucketStatistics::GcBucketStatistics(const GcPailsList &pails, int minAllowedGc, int maxAllowedGc)
    : mMinAllowedGc(minAllowedGc), mMaxAllowedGc(maxAllowedGc), mMeanDepths(101, -1.0)
{
    // 逐字重現 GcBucketStatistics 建構子，**含 off-by-one**：
    // 迴圈第 i 圈取 bucket i+1 的中位數進滑動窗，守門條件檢查 bucket i+1，
    // 但結果寫入 MeanDepths[i]（以 bucket i 為中心的三點平均）。
    double window[3] = { -1.0, -1.0, -1.0 };
    for(int i = 0; i < 100; ++i)              // MeanDepths.length - 1
    {
        const GcPail &pail = pails.buckets()[static_cast<std::size_t>(i + 1)];
        window[0] = window[1];
        window[1] = window[2];
        window[2] = pail.median();
        if(isAllowed(&pail))
        {
            mMeanDepths[static_cast<std::size_t>(i)] = doubleStreamAverage(window, 3);
        }
        else
        {
            mMeanDepths[static_cast<std::size_t>(i)] = -1.0;
        }
    }
    mMeanDepths[100] = -1.0;
}

double GcBucketStatistics::medianReadDepth(const GcPail *pail) const
{
    if(!isAllowed(pail)){ return -1; }
    return mMeanDepths[static_cast<std::size_t>(pail->gc())];
}

}
