#include "Percentile.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cobalt {

namespace {

// KthSelector.select(work, pivots, k)：回傳第 k 小（0-based）。
// 契約與實作細節無關，故用 nth_element。
double select(std::vector<double> &v, std::size_t k)
{
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
    return v[k];
}

}

double percentile(std::vector<double> &values, double quantile)
{
    // EstimationType.evaluate：先檢查範圍
    if(!(quantile > 0.0) || quantile > 100.0)
    {
        throw std::out_of_range("quantile out of range");
    }

    // NaNStrategy.REMOVED：先移除 NaN，n 隨之縮小
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double x){ return std::isnan(x); }),
                 values.end());

    const int n = static_cast<int>(values.size());
    if(n == 0){ return std::nan(""); }               // Java 端在此回 NaN；呼叫端 GCPail.median 先擋掉 n==0

    const double p = quantile / 100.0;

    // EstimationType$1.index（LEGACY）
    double pos;
    if(p == 0.0)      { pos = 0.0; }
    else if(p == 1.0) { pos = static_cast<double>(n); }
    else              { pos = p * (n + 1); }          // 注意 n+1，1-based

    // EstimationType.estimate
    const double fpos = std::floor(pos);
    const int intPos = static_cast<int>(fpos);
    const double dif = pos - fpos;

    if(pos < 1){ return select(values, 0); }                                  // 最小值
    if(pos >= static_cast<double>(n)){ return select(values, static_cast<std::size_t>(n - 1)); }  // 最大值

    const double lower = select(values, static_cast<std::size_t>(intPos - 1));
    const double upper = select(values, static_cast<std::size_t>(intPos));
    // 逐字重現：不可寫成 (lower+upper)/2 或 lower*(1-dif)+upper*dif
    return lower + dif * (upper - lower);
}

}
