#ifndef AMBER_COMMONSMATH_H
#define AMBER_COMMONSMATH_H

// AMBER 用到的 binomial CDF：commons-math3 3.6.1 的
// BinomialDistribution.cumulativeProbability 的等價實作。
//
// 這裡只保留 C++ 標準庫沒有的部分：
//   - 正規化不完全 beta 函數 I_x(a,b)（<cmath> 沒有；std::beta 是**完全** beta，
//     且在此處的 (a,b) 範圍會下溢為 0，不能拿來替代）
//   - 上式所需的連分數展開（Beta.java / ContinuedFraction.java，epsilon 1E-14）
// 其餘一律走標準庫：log / log1p / exp / lgamma。
//
// 對照原始碼：github.com/apache/commons-math tag MATH_3_6_1，
// 以及 amber_v4.3.jar 內 shade 的同版 bytecode。
//
// **不追求與 Java 逐位元相同**：改用標準庫後 CDF 與 Java 有 ulp 級的差異，
// 但在驗證過的樣本上 noiseFloor 與其下游輸出完全不變。

namespace amber { namespace cm3 {

// Beta.regularizedBeta(x, a, b)，epsilon = 1E-14
double regularizedBeta(double x, double a, double b);

// BinomialDistribution(trials, p).cumulativeProbability(x)
double binomialCdf(int numberOfTrials, double probabilityOfSuccess, int x);

}}

#endif
