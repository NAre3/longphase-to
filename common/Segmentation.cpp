#include "Segmentation.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "HumanChromosome.h"

namespace lp {


#include "Centromeres38.inc"

int centromereOf(const std::string &chromosome)
{
    const std::string stripped = stripChrPrefix(chromosome);
    for(const auto &entry : CENTROMERES_38){
        if(stripped == entry.first){
            return entry.second;
        }
    }
    throw std::runtime_error("no centromere for chromosome " + chromosome);
}

// HumanChromosome 的宣告順序（HumanChromosome.java:19-42）：_1.._22, _X, _Y。
// ChrArm.compareTo 先比染色體的 enum 序、再比 arm（P < Q），SegmentsFile 與 CP-A8/A9
// 的輸出順序皆依此排序。
int chromosomeOrdinal(const std::string &shortName)
{
    if(shortName == "X"){ return 22; }
    if(shortName == "Y"){ return 23; }
    return std::stoi(shortName) - 1;
}

// Doubles.median（Doubles.java:90-104）：排序後取中位數，偶數長度取中間兩數平均
static double median(std::vector<double> values)
{
    if(values.empty()){
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if(n % 2 == 0){
        return (values[n / 2] + values[n / 2 - 1]) / 2;
    }
    return values[n / 2];
}

// Doubles.mean（Doubles.java:75-88）：依序累加後除以長度。
// 累加順序影響最後一位。【2026-09-06 更正理由】原註寫「必須照抄順序」，那是追求與 Java
// 逐位元相同時的說法；現行理由是本程式自己的輸出要穩定——這個值會經 DecimalFormat("#.####")
// 寫進 .pcf 的 MeanRatio，由 PURPLE 讀取。此迴圈與 std::accumulate 等價，改寫無益亦無害；
// 不要改成 pairwise 或 Kahan 求和，那會讓輸出隨實作細節變動。
double mean(const double *values, std::size_t count)
{
    if(count == 0){
        return std::numeric_limits<double>::quiet_NaN();
    }
    double sum = 0.0;
    for(std::size_t i = 0; i < count; ++i){
        sum += values[i];
    }
    return sum / count;
}

// Stats.medianAbsoluteDeviation（Stats.java）
static double medianAbsoluteDeviation(const std::vector<double> &values)
{
    const double med = median(values);
    std::vector<double> deviations(values.size());
    for(std::size_t i = 0; i < values.size(); ++i){
        deviations[i] = std::fabs(values[i] - med);
    }
    return median(std::move(deviations)) * 1.4826;
}

// Runmed.med3
static double med3(double a, double b, double c)
{
    if(a < b){
        if(c < b){
            return std::fmax(a, c);
        }
        return b;
    }
    if(c < a){
        return std::fmax(c, b);
    }
    return a;
}

// Runmed.medianOdd：排序後取第 (n+1)/2 小
static double medianOdd(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[(values.size() + 1) / 2 - 1];
}

// WindowedMedian（WindowedMedian.java）。
// Java 以雙 heap 維護視窗，getMedian() 回傳 maxHeap.peek()，即視窗內第
// ceil(w/2) 小的值。視窗大小為奇數且不超過資料長度時，該值即真正的中位數。
// 此處以「第 k 小」的形式實作，與 heap 的機制無關但輸出相同。
static std::vector<double> windowedMedians(const std::vector<double> &data, int windowSize)
{
    const int n = static_cast<int>(data.size());
    const int halfWindow = windowSize / 2;
    std::vector<double> result(data.size(), 0.0);

    if(n == 0){
        return result;
    }

    std::vector<double> window;
    window.reserve(windowSize);

    const int initial = std::min(windowSize, n);
    for(int i = 0; i < initial; ++i){
        window.push_back(data[i]);
    }

    const auto kthSmallest = [](std::vector<double> values, int k){
        std::nth_element(values.begin(), values.begin() + (k - 1), values.end());
        return values[k - 1];
    };

    if(halfWindow < n){
        // maxHeap.size() = ceil(size/2)
        const int k = (static_cast<int>(window.size()) + 1) / 2;
        result[halfWindow] = kthSmallest(window, k);
    }

    for(int i = windowSize; i < n; ++i){
        // 移除 data[i - windowSize]、加入 data[i]
        const double leaving = data[i - windowSize];
        const auto it = std::find(window.begin(), window.end(), leaving);
        if(it != window.end()){
            window.erase(it);
        }
        window.push_back(data[i]);

        const int k = (static_cast<int>(window.size()) + 1) / 2;
        result[i - halfWindow] = kthSmallest(window, k);
    }

    // 前後 halfWindow 個位置直接複製原值
    for(int i = 0; i < halfWindow && i < n; ++i){
        result[i] = data[i];
    }
    for(int i = n - halfWindow; i < n; ++i){
        if(i >= 0){
            result[i] = data[i];
        }
    }

    return result;
}

// Runmed.smoothEnds
static std::vector<double> smoothEnds(const std::vector<double> &y, int k)
{
    const int halfK = k / 2;
    std::vector<double> sm = y;

    if(halfK < 1){
        return sm;
    }

    const int n = static_cast<int>(y.size());

    if(halfK >= 2){
        sm[1] = med3(y[0], y[1], y[2]);
        sm[n - 2] = med3(y[n - 1], y[n - 2], y[n - 3]);

        if(halfK >= 3){
            for(int i = 3; i <= halfK; ++i){
                const int j = 2 * i - 1;
                sm[i - 1] = medianOdd(std::vector<double>(y.begin(), y.begin() + j));
                sm[n - i] = medianOdd(std::vector<double>(y.begin() + (n - j), y.begin() + n));
            }
        }
    }

    sm[0] = med3(y[0], sm[1], 3 * sm[1] - 2 * sm[2]);
    sm[n - 1] = med3(y[n - 1], sm[n - 2], 3 * sm[n - 2] - 2 * sm[n - 3]);

    return sm;
}

// Gamma.filterWidth
static int filterWidth(int n)
{
    return n >= 51 ? 51 : (n % 2 == 0 ? n - 1 : n);
}


std::vector<double> runmed(const std::vector<double> &data, int k, bool smooth)
{
    if(k < 1 || k % 2 != 1){
        throw std::runtime_error("runmed window must be odd and >= 1");
    }
    if(data.size() < 2){
        return data;
    }
    const std::vector<double> medians = windowedMedians(data, k);
    return smooth ? smoothEnds(medians, k) : medians;
}

// Gamma（Gamma.java）
double gammaSegmentPenalty(const std::vector<double> &y, double gamma, bool normalise)
{
    if(y.empty()){
        throw std::runtime_error("Input array must not be empty");
    }

    if(!normalise){
        return gamma;
    }

    const std::vector<double> runningMedians = runmed(y, filterWidth(static_cast<int>(y.size())), true);

    std::vector<double> diffs(y.size());
    for(std::size_t i = 0; i < y.size(); ++i){
        diffs[i] = y[i] - runningMedians[i];
    }

    const double sd = medianAbsoluteDeviation(diffs);
    return sd == 0 ? 0.01 * gamma : sd * sd * gamma;
}

// Segmenter（Segmenter.java）：最小成本分段的動態規劃。C++ 標準庫沒有對應設施，
// 這是演算法本身。比較方式（嚴格 <）與 Double.MAX_VALUE 的初值決定成本相同時選哪個切點，
// 改成 <= 會選到不同切點——那是**分段結果不同**，不是最後一位不同，與浮點精度無關。
Fit segment(const std::vector<double> &y, double segmentPenalty)
{
    const std::size_t n = y.size();

    std::vector<double> cumulativeSums(n);
    std::vector<double> cumulativeSquaredSums(n);
    if(n > 0){
        cumulativeSums[0] = y[0];
        cumulativeSquaredSums[0] = y[0] * y[0];
        for(std::size_t i = 1; i < n; ++i){
            cumulativeSums[i] = cumulativeSums[i - 1] + y[i];
            cumulativeSquaredSums[i] = cumulativeSquaredSums[i - 1] + (y[i] * y[i]);
        }
    }

    std::vector<double> leastCostEndingJustBefore(n + 1, 0.0);
    std::vector<int> leastCostSegmentEndpoints(n, 0);

    for(std::size_t end = 0; end < n; ++end){
        double minCost = std::numeric_limits<double>::max(); // Double.MAX_VALUE
        int endOfPreviousSegmentForLeastCost = 0;

        for(std::size_t start = 0; start <= end; ++start){
            const double sums = start == 0
                    ? cumulativeSums[end]
                    : cumulativeSums[end] - cumulativeSums[start - 1];
            const double sumSquared = sums * sums;
            const double sumOfSquares = start == 0
                    ? cumulativeSquaredSums[end]
                    : cumulativeSquaredSums[end] - cumulativeSquaredSums[start - 1];

            const int segmentLength = static_cast<int>(end - start + 1);
            const double segmentCost = sumOfSquares - (sumSquared / segmentLength);
            const double cost = leastCostEndingJustBefore[start] + segmentPenalty + segmentCost;

            if(cost < minCost){
                minCost = cost;
                endOfPreviousSegmentForLeastCost = static_cast<int>(start) - 1;
            }
        }

        leastCostEndingJustBefore[end + 1] = minCost;
        leastCostSegmentEndpoints[end] = endOfPreviousSegmentForLeastCost;
    }

    std::vector<int> segmentEndpoints;
    int lastSegmentEndpoint = static_cast<int>(n) - 1;
    while(lastSegmentEndpoint >= 0){
        segmentEndpoints.push_back(lastSegmentEndpoint);
        lastSegmentEndpoint = leastCostSegmentEndpoints[lastSegmentEndpoint];
    }
    std::reverse(segmentEndpoints.begin(), segmentEndpoints.end());

    Fit fit;
    int start = 0;
    for(int endpoint : segmentEndpoints){
        fit.lengths.push_back(endpoint - start + 1);
        fit.startPositions.push_back(start);
        start = endpoint + 1;
    }

    return fit;
}

}
