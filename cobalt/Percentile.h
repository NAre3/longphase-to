#ifndef COBALT_PERCENTILE_H
#define COBALT_PERCENTILE_H

#include <vector>

namespace cobalt {

// commons-math3 DescriptiveStatistics.getPercentile(p) 的移植。
// 組態由 Percentile(double) 建構子硬寫：EstimationType.LEGACY + NaNStrategy.REMOVED
// + KthSelector(MedianOf3PivotingStrategy)。完整依據見 behaviour-contract-percentile.md。
//
// 內插規則屬演算法邏輯，逐字重現；select 的第 k 小以 std::nth_element 達成同一契約
// （該契約與 introselect 的樞紐細節無關）。
//
// values 會被就地重排（與 Java 端在 work array 上操作相同，不影響呼叫端語義）。
double percentile(std::vector<double> &values, double quantile);

}

#endif
