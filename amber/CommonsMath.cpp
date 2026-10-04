#include "CommonsMath.h"

#include <cmath>
#include <limits>
#include <stdexcept>

// 本檔對應 commons-math3 3.6.1 的 Beta.java 與 ContinuedFraction.java。
//
// 2026-09-06：原本連 FastMath 與 Gamma 的高精度 logBeta 機制一併照抄，
// 目的是與 Java 逐位元相同。既然邏輯等價即可，改為：
//   - log / log1p / exp  → 標準庫
//   - logBeta(p,q)       → lgamma(p) + lgamma(q) - lgamma(p+q)
//     （commons-math 的 deltaMinusDeltaSum / invGamma1pm1 / logGamma1p 機制整組移除）
// 連分數展開保留：正規化不完全 beta 標準庫沒有，那是演算法本身，不是重造輪子。

namespace amber { namespace cm3 {

namespace {

// Beta.java: DEFAULT_EPSILON
constexpr double DEFAULT_EPSILON = 1E-14;

// Precision.equals(x, 0.0, small) 在此僅用於與 0 的比較，等價於 |x| <= small
inline bool nearZero(double value, double small)
{
    return std::fabs(value) <= small;
}

// Beta.logBeta(p, q)。commons-math 為了大參數下的精度用了一整套 delta 修正；
// 此處採用等價的 lgamma 形式，差異為 ulp 級（實測見 math_provenance/README.md）。
double logBeta(double p, double q)
{
    if(std::isnan(p) || std::isnan(q) || p <= 0.0 || q <= 0.0){
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::lgamma(p) + std::lgamma(q) - std::lgamma(p + q);
}

// Beta.java 內的 ContinuedFraction 求值（getA 恆為 1，getB 依奇偶分支）
double betaFractionEvaluate(double x, double a, double b, double epsilon)
{
    const double small = 1e-50;

    double hPrev = 1.0; // getA(0, x)
    if(nearZero(hPrev, small)){
        hPrev = small;
    }

    int n = 1;
    double dPrev = 0.0;
    double cPrev = hPrev;
    double hN = hPrev;

    while(true){
        const double fa = 1.0; // getA(n, x)

        double fb;
        if(n % 2 == 0){
            const double m = n / 2.0;
            fb = (m * (b - m) * x) / ((a + (2 * m) - 1) * (a + (2 * m)));
        }else{
            const double m = (n - 1.0) / 2.0;
            fb = -((a + m) * (a + b + m) * x) / ((a + (2 * m)) * (a + (2 * m) + 1.0));
        }

        double dN = fa + fb * dPrev;
        if(nearZero(dN, small)){
            dN = small;
        }
        double cN = fa + fb / cPrev;
        if(nearZero(cN, small)){
            cN = small;
        }

        dN = 1 / dN;
        const double deltaN = cN * dN;
        hN = hPrev * deltaN;

        if(std::isinf(hN) || std::isnan(hN)){
            throw std::runtime_error("continued fraction diverged");
        }

        if(std::fabs(deltaN - 1.0) < epsilon){
            break;
        }

        dPrev = dN;
        cPrev = cN;
        hPrev = hN;
        ++n;
    }

    return hN;
}

}

double regularizedBeta(double x, double a, double b)
{
    if(std::isnan(x) || std::isnan(a) || std::isnan(b) || x < 0 || x > 1 || a <= 0 || b <= 0){
        return std::numeric_limits<double>::quiet_NaN();
    }

    if(x > (a + 1) / (2 + b + a) && 1 - x <= (b + 1) / (2 + b + a)){
        return 1 - regularizedBeta(1 - x, b, a);
    }

    return std::exp((a * std::log(x)) + (b * std::log1p(-x)) - std::log(a) - logBeta(a, b))
            * 1.0 / betaFractionEvaluate(x, a, b, DEFAULT_EPSILON);
}

double binomialCdf(int numberOfTrials, double probabilityOfSuccess, int x)
{
    if(x < 0){
        return 0.0;
    }
    if(x >= numberOfTrials){
        return 1.0;
    }
    return 1.0 - regularizedBeta(probabilityOfSuccess, x + 1.0, numberOfTrials - x);
}

}}
