#include "Segmentation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>

#include "CpDump.h"
#include "HumanChromosome.h"

namespace amber {

namespace {

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
double median(std::vector<double> values)
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
// 累加順序影響最後一位，故必須照抄順序，不可用其他求和方式。
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
double medianAbsoluteDeviation(const std::vector<double> &values)
{
    const double med = median(values);
    std::vector<double> deviations(values.size());
    for(std::size_t i = 0; i < values.size(); ++i){
        deviations[i] = std::fabs(values[i] - med);
    }
    return median(std::move(deviations)) * 1.4826;
}

// Runmed.med3
double med3(double a, double b, double c)
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
double medianOdd(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[(values.size() + 1) / 2 - 1];
}

// WindowedMedian（WindowedMedian.java）。
// Java 以雙 heap 維護視窗，getMedian() 回傳 maxHeap.peek()，即視窗內第
// ceil(w/2) 小的值。視窗大小為奇數且不超過資料長度時，該值即真正的中位數。
// 此處以「第 k 小」的形式實作，與 heap 的機制無關但輸出相同。
std::vector<double> windowedMedians(const std::vector<double> &data, int windowSize)
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
std::vector<double> smoothEnds(const std::vector<double> &y, int k)
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
int filterWidth(int n)
{
    return n >= 51 ? 51 : (n % 2 == 0 ? n - 1 : n);
}

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

namespace {

struct Fit
{
    std::vector<int> lengths;
    std::vector<int> startPositions;
};

// Segmenter（Segmenter.java）：最小成本分段的動態規劃。
// 累加順序、比較方式（嚴格 <）與 Double.MAX_VALUE 的初值都照抄；
// 這些決定了在成本相同時選哪一個切點。
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

SegmentationResult segmentBafs(const std::vector<AmberBAF> &bafs)
{
    SegmentationResult result;

    // PerArmSegmenter 建構（PerArmSegmenter.java:39-52）：
    // 走訪各染色體的 BAF，value(baf) = tumorModifiedBAF，>= 0.0 者依 (染色體, arm) 分組。
    // tumorModifiedBAF = 0.5 + |baf - 0.5| >= 0.5，故實務上全部入選。
    std::map<std::string, std::size_t> armIndex;

    for(const AmberBAF &baf : bafs){
        const double value = baf.tumorModifiedBAF();
        if(!(value >= 0.0)){
            continue;
        }

        const std::string shortName = stripChrPrefix(baf.chromosome);
        const char arm = baf.position < centromereOf(baf.chromosome) ? 'P' : 'Q';
        const std::string armId = shortName + "_" + arm;

        auto found = armIndex.find(armId);
        if(found == armIndex.end()){
            armIndex.emplace(armId, result.arms.size());
            ArmSegments segments;
            segments.armId = armId;
            segments.chromosomeShort = shortName;
            segments.chromosomeRank = chromosomeOrdinal(shortName);
            segments.arm = arm;
            result.arms.push_back(std::move(segments));
            found = armIndex.find(armId);
        }

        result.arms[found->second].values.push_back(value);
    }

    // ChrArm.compareTo：先比染色體 enum 序，再比 arm（P < Q）
    std::sort(result.arms.begin(), result.arms.end(), [](const ArmSegments &a, const ArmSegments &b){
        if(a.chromosomeRank != b.chromosomeRank){
            return a.chromosomeRank < b.chromosomeRank;
        }
        return a.arm < b.arm;
    });

    int totalCount = 0;
    for(const ArmSegments &arm : result.arms){
        totalCount += static_cast<int>(arm.values.size());
    }
    result.totalCount = totalCount;
    result.penaltyMode = totalCount < UNIFORM_PENALTY_THRESHOLD ? "uniform" : "per-arm-gamma";

    if(result.penaltyMode == "uniform"){
        // PerArmSegmenter.java:55-72 的 uniform 分支。其 allRatios 的組裝順序取自
        // HashMap.keySet()，跨語言不保證。本研究的資料不走此分支（spec U5 已關閉），
        // 故不實作——寧可明確失敗，也不要用一個未經驗證的順序默默算出結果。
        throw std::runtime_error(
                "uniform penalty branch is not implemented: its allRatios assembly order comes from "
                "HashMap.keySet() and is not reproducible across languages (see spec U5). "
                "totalCount=" + std::to_string(totalCount));
    }

    // 需要每個 arm 內的 BAF（依原順序）以取得 segment 的起訖座標
    std::map<std::string, std::vector<const AmberBAF *>> bafsByArm;
    for(const AmberBAF &baf : bafs){
        if(!(baf.tumorModifiedBAF() >= 0.0)){
            continue;
        }
        const std::string shortName = stripChrPrefix(baf.chromosome);
        const char arm = baf.position < centromereOf(baf.chromosome) ? 'P' : 'Q';
        bafsByArm[shortName + "_" + arm].push_back(&baf);
    }

    for(ArmSegments &arm : result.arms){
        if(arm.values.empty()){
            continue;
        }

        // Segmenter 以 GammaPenaltyCalculator(gamma, true) 逐 arm 計算 penalty
        const double penalty = gammaSegmentPenalty(arm.values, BAF_SEGMENTATION_GAMMA, true);
        const Fit fit = segment(arm.values, penalty);

        // ChromosomeArmSegments 建構（ChromosomeArmSegments.java:14-30）：
        // MeanRatio 取 rawValues 該段的平均（Doubles.mean），**不是** PiecewiseConstantFit
        // 裡那個已四捨五入到三位小數的 means。
        const std::vector<const AmberBAF *> &armBafs = bafsByArm[arm.armId];
        std::size_t ratiosIndex = 0;

        for(std::size_t i = 0; i < fit.lengths.size(); ++i){
            const int count = fit.lengths[i];
            const double meanRatio = mean(arm.values.data() + ratiosIndex, static_cast<std::size_t>(count));

            const AmberBAF *startRatio = armBafs[ratiosIndex];
            ratiosIndex += static_cast<std::size_t>(count);
            const AmberBAF *endRatio = armBafs[ratiosIndex - 1];

            PcfSegment segment;
            segment.chromosome = arm.chromosomeShort;
            segment.start = startRatio->position;
            segment.end = endRatio->position; // AbsoluteSegments.segmentEnd = endRatio.position()
            segment.meanRatio = meanRatio;
            arm.segments.push_back(std::move(segment));
        }
    }

    return result;
}

namespace {

// 對應 CobaltRatioFile.FORMAT = DecimalFormat("#.####")：
// 最多四位小數、去掉尾端的 0 與孤立的小數點。
std::string format4Significant(double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.4f", value);
    std::string s = buffer;

    if(s.find('.') != std::string::npos){
        while(!s.empty() && s.back() == '0'){
            s.pop_back();
        }
        if(!s.empty() && s.back() == '.'){
            s.pop_back();
        }
    }

    return s;
}

}

void writeSegmentsFile(const std::string &path, const SegmentationResult &result)
{
    std::ofstream out(path);
    if(!out){
        throw std::runtime_error("unable to open for writing: " + path);
    }

    out << "Chromosome\tStart\tEnd\tMeanRatio\n";

    for(const ArmSegments &arm : result.arms){
        for(const PcfSegment &segment : arm.segments){
            // SegmentsFile 以 genomeVersion.versionedChromosome 加上 chr 前綴
            out << "chr" << arm.chromosomeShort << '\t'
                << segment.start << '\t'
                << segment.end << '\t'
                << format4Significant(segment.meanRatio) << '\n';
        }
    }
}


void writeSegmentationCheckpoints(const SegmentationResult &result)
{
    if(!CpDump::enabled()){
        return;
    }

    // CP-A8：分段輸入與 penalty（對應 patch_segmenter.py 插在 PerArmSegmenter 建構子尾端）
    std::vector<CpDump::Row> a8;
    a8.push_back({"", 0, "totalCount\t" + std::to_string(result.totalCount)});
    a8.push_back({"", 0, "uniformPenaltyThreshold\t" + std::to_string(result.uniformPenaltyThreshold)});
    a8.push_back({"", 0, "penaltyMode\t" + result.penaltyMode});
    a8.push_back({"", 0, "gamma\t" + CpDump::num(result.gamma)});

    for(const ArmSegments &arm : result.arms){
        a8.push_back({"", 0, "armCount." + arm.armId + "\t" + std::to_string(arm.values.size())});
        for(std::size_t i = 0; i < arm.values.size(); ++i){
            // BAFSegmenter 的 valuesForSegmentation 與 rawValues 相同（BAFSegmenter.java:47-58）
            a8.push_back({"", 0, "armValue." + arm.armId + "." + std::to_string(i) + "\t"
                    + CpDump::num(arm.values[i]) + "\t" + CpDump::num(arm.values[i])});
        }
    }
    CpDump::write("CP-A8", "field\tvalue\trawValue", a8);

    // CP-A9：分段結果
    std::vector<CpDump::Row> a9;
    for(const ArmSegments &arm : result.arms){
        int idx = 0;
        for(const PcfSegment &segment : arm.segments){
            a9.push_back({"", 0, arm.armId + "\t" + std::to_string(idx++) + "\t"
                    + segment.chromosome + "\t" + std::to_string(segment.start) + "\t"
                    + std::to_string(segment.end) + "\t" + CpDump::num(segment.meanRatio)});
        }
    }
    CpDump::write("CP-A9", "arm\tidx\tchromosome\tstart\tend\tmeanRatio", a9);
}

}
