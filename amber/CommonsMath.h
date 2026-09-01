#ifndef AMBER_COMMONSMATH_H
#define AMBER_COMMONSMATH_H

// Apache Commons Math 3.6.1 的 BinomialDistribution.cumulativeProbability 在 C++ 端的等價實作。
//
// 為什麼要照抄到這個程度：AMBER 的 peak 捕捉判定是拿 CDF 去比 0.16 / 0.84 兩個門檻
// （CandidatePeak.java:88-93）。門檻判定是布林值，CDF 的最後一個位元不同就可能讓某個點
// 一邊算進去、一邊沒算進去。因此連續分數的收斂條件（epsilon 1E-14）與 logBeta 的
// 分支結構都屬於契約的一部分，不只是「算出同一個數學函數」。
//
// commons-math3 3.6.1 已 shade 在 amber_v4.3.jar 內
// （org/apache/commons/math3/{distribution/BinomialDistribution,special/Beta,special/Gamma,
//   util/ContinuedFraction}.class，時間戳 2016-03-17），該 bytecode 即產生參考值的實體。
//
// 對照原始碼：github.com/apache/commons-math tag MATH_3_6_1

namespace amber { namespace cm3 {

// Beta.regularizedBeta(x, a, b) —— 使用 DEFAULT_EPSILON = 1E-14、maxIterations = INT_MAX
double regularizedBeta(double x, double a, double b);

// Beta.logBeta(p, q)
double logBeta(double p, double q);

// Gamma.logGamma(x)
double logGamma(double x);

// Gamma.logGamma1p(x)，-0.5 <= x <= 1.5
double logGamma1p(double x);

// Gamma.invGamma1pm1(x)，-0.5 <= x <= 1.5
double invGamma1pm1(double x);

// BinomialDistribution(trials, p).cumulativeProbability(x)
double binomialCdf(int numberOfTrials, double probabilityOfSuccess, int x);

}}

#endif
