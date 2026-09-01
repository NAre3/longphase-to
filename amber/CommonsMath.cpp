#include "CommonsMath.h"

#include "FastMath.h"

#include <cmath>
#include <limits>
#include <stdexcept>

// 本檔逐行對應 commons-math3 3.6.1 的 Beta.java / Gamma.java / ContinuedFraction.java。
// FastMath.log / log1p / exp 一律走移植版 amber::fm（見 FastMath.h）。
// 先前以 std:: 版本實作時，61366 組對照中有 215 組與 Java 不逐位元相同（1~32 ulp），
// 因為 commons-math 的 FastMath 與系統 libm 本來就不同——實測 FastMath 與 java.lang.Math
// 在 20 萬組取樣中 log1p 差 13142 組、exp 差 457 組、log 差 21 組。
// 對照工具：tools/CdfConformance.java（Java 端基準）與 tools/cdf_conformance.cpp（比對）。

namespace amber { namespace cm3 {

namespace {

// Gamma.java 的 INV_GAMMA1P_M1_* 常數，由原始碼機械抽出（未經人工轉錄）
constexpr double INV_GAMMA1P_M1_A0 = .611609510448141581788E-08;
constexpr double INV_GAMMA1P_M1_A1 = .624730830116465516210E-08;
constexpr double INV_GAMMA1P_M1_B1 = .203610414066806987300E+00;
constexpr double INV_GAMMA1P_M1_B2 = .266205348428949217746E-01;
constexpr double INV_GAMMA1P_M1_B3 = .493944979382446875238E-03;
constexpr double INV_GAMMA1P_M1_B4 = -.851419432440314906588E-05;
constexpr double INV_GAMMA1P_M1_B5 = -.643045481779353022248E-05;
constexpr double INV_GAMMA1P_M1_B6 = .992641840672773722196E-06;
constexpr double INV_GAMMA1P_M1_B7 = -.607761895722825260739E-07;
constexpr double INV_GAMMA1P_M1_B8 = .195755836614639731882E-09;
constexpr double INV_GAMMA1P_M1_P0 = .6116095104481415817861E-08;
constexpr double INV_GAMMA1P_M1_P1 = .6871674113067198736152E-08;
constexpr double INV_GAMMA1P_M1_P2 = .6820161668496170657918E-09;
constexpr double INV_GAMMA1P_M1_P3 = .4686843322948848031080E-10;
constexpr double INV_GAMMA1P_M1_P4 = .1572833027710446286995E-11;
constexpr double INV_GAMMA1P_M1_P5 = -.1249441572276366213222E-12;
constexpr double INV_GAMMA1P_M1_P6 = .4343529937408594255178E-14;
constexpr double INV_GAMMA1P_M1_Q1 = .3056961078365221025009E+00;
constexpr double INV_GAMMA1P_M1_Q2 = .5464213086042296536016E-01;
constexpr double INV_GAMMA1P_M1_Q3 = .4956830093825887312020E-02;
constexpr double INV_GAMMA1P_M1_Q4 = .2692369466186361192876E-03;
constexpr double INV_GAMMA1P_M1_C = -.422784335098467139393487909917598E+00;
constexpr double INV_GAMMA1P_M1_C0 = .577215664901532860606512090082402E+00;
constexpr double INV_GAMMA1P_M1_C1 = -.655878071520253881077019515145390E+00;
constexpr double INV_GAMMA1P_M1_C2 = -.420026350340952355290039348754298E-01;
constexpr double INV_GAMMA1P_M1_C3 = .166538611382291489501700795102105E+00;
constexpr double INV_GAMMA1P_M1_C4 = -.421977345555443367482083012891874E-01;
constexpr double INV_GAMMA1P_M1_C5 = -.962197152787697356211492167234820E-02;
constexpr double INV_GAMMA1P_M1_C6 = .721894324666309954239501034044657E-02;
constexpr double INV_GAMMA1P_M1_C7 = -.116516759185906511211397108401839E-02;
constexpr double INV_GAMMA1P_M1_C8 = -.215241674114950972815729963053648E-03;
constexpr double INV_GAMMA1P_M1_C9 = .128050282388116186153198626328164E-03;
constexpr double INV_GAMMA1P_M1_C10 = -.201348547807882386556893914210218E-04;
constexpr double INV_GAMMA1P_M1_C11 = -.125049348214267065734535947383309E-05;
constexpr double INV_GAMMA1P_M1_C12 = .113302723198169588237412962033074E-05;
constexpr double INV_GAMMA1P_M1_C13 = -.205633841697760710345015413002057E-06;

// 註：Gamma.java 的 HALF_LOG_2_PI 只用於 logGamma 的 Lanczos 分支，該分支在此不可達
// （見 logGamma 末端的 throw），故不在此定義，避免看起來像已驗證的常數。

// Beta.java: HALF_LOG_TWO_PI
constexpr double HALF_LOG_TWO_PI = 0.9189385332046727;

// Beta.java: DELTA
constexpr double DELTA[] = {
    .833333333333333333333333333333E-01,
    -.277777777777777777777777752282E-04,
    .793650793650793650791732130419E-07,
    -.595238095238095232389839236182E-09,
    .841750841750832853294451671990E-11,
    -.191752691751854612334149171243E-12,
    .641025640510325475730918472625E-14,
    -.295506514125338232839867823991E-15,
    .179643716359402238723287696452E-16,
    -.139228964661627791231203060395E-17,
    .133802855014020915603275339093E-18,
    -.154246009867966094273710216533E-19,
    .197701992980957427278370133333E-20,
    -.234065664793997056856992426667E-21,
    .171348014966398575409015466667E-22
};
constexpr int DELTA_LEN = sizeof(DELTA) / sizeof(DELTA[0]);

// Beta.java: DEFAULT_EPSILON
constexpr double DEFAULT_EPSILON = 1E-14;

// Precision.equals(x, 0.0, small) 在此僅用於與 0 的比較，等價於 |x| <= small
inline bool nearZero(double value, double small)
{
    return std::fabs(value) <= small;
}

// Beta.java:deltaMinusDeltaSum（要求 0 <= a <= b 且 b >= 10）
double deltaMinusDeltaSum(double a, double b)
{
    if(a < 0 || a > b || b < 10){
        throw std::runtime_error("deltaMinusDeltaSum out of range");
    }

    const double h = a / b;
    const double p = h / (1.0 + h);
    const double q = 1.0 / (1.0 + h);
    const double q2 = q * q;

    double s[DELTA_LEN];
    s[0] = 1.0;
    for(int i = 1; i < DELTA_LEN; ++i){
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }

    const double sqrtT = 10.0 / b;
    const double t = sqrtT * sqrtT;
    double w = DELTA[DELTA_LEN - 1] * s[DELTA_LEN - 1];
    for(int i = DELTA_LEN - 2; i >= 0; --i){
        w = t * w + DELTA[i] * s[i];
    }

    return w * p / b;
}

// Beta.java:sumDeltaMinusDeltaSum（要求 p >= 10 且 q >= 10）
double sumDeltaMinusDeltaSum(double p, double q)
{
    if(p < 10.0 || q < 10.0){
        throw std::runtime_error("sumDeltaMinusDeltaSum out of range");
    }

    const double a = std::fmin(p, q);
    const double b = std::fmax(p, q);
    const double sqrtT = 10.0 / a;
    const double t = sqrtT * sqrtT;

    double z = DELTA[DELTA_LEN - 1];
    for(int i = DELTA_LEN - 2; i >= 0; --i){
        z = t * z + DELTA[i];
    }

    return z / a + deltaMinusDeltaSum(a, b);
}

// Beta.java:logGammaSum（要求 1 <= a,b <= 2）
double logGammaSum(double a, double b)
{
    if(a < 1.0 || a > 2.0 || b < 1.0 || b > 2.0){
        throw std::runtime_error("logGammaSum out of range");
    }

    const double x = (a - 1.0) + (b - 1.0);

    if(x <= 0.5){
        return logGamma1p(1.0 + x);
    }
    if(x <= 1.5){
        return logGamma1p(x) + fm::log1p(x);
    }
    return logGamma1p(x - 1.0) + fm::log(x * (1.0 + x));
}

// Beta.java:logGammaMinusLogGammaSum（要求 a >= 0 且 b >= 10）
double logGammaMinusLogGammaSum(double a, double b)
{
    if(a < 0.0 || b < 10.0){
        throw std::runtime_error("logGammaMinusLogGammaSum out of range");
    }

    double d;
    double w;
    if(a <= b){
        d = b + (a - 0.5);
        w = deltaMinusDeltaSum(a, b);
    }else{
        d = a + (b - 0.5);
        w = deltaMinusDeltaSum(b, a);
    }

    const double u = d * fm::log1p(a / b);
    const double v = a * (fm::log(b) - 1.0);

    return u <= v ? (w - u) - v : (w - v) - u;
}

// ContinuedFraction.evaluate 針對 regularizedBeta 的 getA/getB 特化
// （Beta.java:198-222 的匿名子類；getA 恆為 1.0）
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

double invGamma1pm1(double x)
{
    if(x < -0.5 || x > 1.5){
        throw std::runtime_error("invGamma1pm1 out of range");
    }

    double ret;
    const double t = x <= 0.5 ? x : (x - 0.5) - 0.5;

    if(t < 0.0){
        const double a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
        double b = INV_GAMMA1P_M1_B8;
        b = INV_GAMMA1P_M1_B7 + t * b;
        b = INV_GAMMA1P_M1_B6 + t * b;
        b = INV_GAMMA1P_M1_B5 + t * b;
        b = INV_GAMMA1P_M1_B4 + t * b;
        b = INV_GAMMA1P_M1_B3 + t * b;
        b = INV_GAMMA1P_M1_B2 + t * b;
        b = INV_GAMMA1P_M1_B1 + t * b;
        b = 1.0 + t * b;

        double c = INV_GAMMA1P_M1_C13 + t * (a / b);
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C + t * c;

        if(x > 0.5){
            ret = t * c / x;
        }else{
            ret = x * ((c + 0.5) + 0.5);
        }
    }else{
        double p = INV_GAMMA1P_M1_P6;
        p = INV_GAMMA1P_M1_P5 + t * p;
        p = INV_GAMMA1P_M1_P4 + t * p;
        p = INV_GAMMA1P_M1_P3 + t * p;
        p = INV_GAMMA1P_M1_P2 + t * p;
        p = INV_GAMMA1P_M1_P1 + t * p;
        p = INV_GAMMA1P_M1_P0 + t * p;

        double q = INV_GAMMA1P_M1_Q4;
        q = INV_GAMMA1P_M1_Q3 + t * q;
        q = INV_GAMMA1P_M1_Q2 + t * q;
        q = INV_GAMMA1P_M1_Q1 + t * q;
        q = 1.0 + t * q;

        double c = INV_GAMMA1P_M1_C13 + (p / q) * t;
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C0 + t * c;

        if(x > 0.5){
            ret = (t / x) * ((c - 0.5) - 0.5);
        }else{
            ret = x * c;
        }
    }

    return ret;
}

double logGamma1p(double x)
{
    if(x < -0.5 || x > 1.5){
        throw std::runtime_error("logGamma1p out of range");
    }
    return -fm::log1p(invGamma1pm1(x));
}

double logGamma(double x)
{
    if(std::isnan(x) || x <= 0.0){
        return std::numeric_limits<double>::quiet_NaN();
    }

    if(x < 0.5){
        return logGamma1p(x) - fm::log(x);
    }
    if(x <= 2.5){
        return logGamma1p((x - 0.5) - 0.5);
    }
    if(x <= 8.0){
        const int n = static_cast<int>(std::floor(x - 1.5));
        double prod = 1.0;
        for(int i = 1; i <= n; ++i){
            prod *= x - i;
        }
        return logGamma1p(x - (n + 1)) + fm::log(prod);
    }

    // Lanczos 分支。本研究的呼叫皆為 logBeta 內的 Gamma.logGamma(ared)，ared 落在 (1, 2.5]，
    // 故此分支不應被觸發；若觸發代表前提錯誤，寧可明確失敗也不要靜默用近似值。
    throw std::runtime_error("logGamma lanczos branch reached (unexpected for this study)");
}

double logBeta(double p, double q)
{
    if(std::isnan(p) || std::isnan(q) || p <= 0.0 || q <= 0.0){
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double a = std::fmin(p, q);
    const double b = std::fmax(p, q);

    if(a >= 10.0){
        const double w = sumDeltaMinusDeltaSum(a, b);
        const double h = a / b;
        const double c = h / (1.0 + h);
        const double u = -(a - 0.5) * fm::log(c);
        const double v = b * fm::log1p(h);
        if(u <= v){
            return (((-0.5 * fm::log(b) + HALF_LOG_TWO_PI) + w) - u) - v;
        }
        return (((-0.5 * fm::log(b) + HALF_LOG_TWO_PI) + w) - v) - u;
    }

    if(a > 2.0){
        if(b > 1000.0){
            const int n = static_cast<int>(std::floor(a - 1.0));
            double prod = 1.0;
            double ared = a;
            for(int i = 0; i < n; ++i){
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            return (fm::log(prod) - n * fm::log(b))
                    + (logGamma(ared) + logGammaMinusLogGammaSum(ared, b));
        }

        double prod1 = 1.0;
        double ared = a;
        while(ared > 2.0){
            ared -= 1.0;
            const double h = ared / b;
            prod1 *= h / (1.0 + h);
        }

        if(b < 10.0){
            double prod2 = 1.0;
            double bred = b;
            while(bred > 2.0){
                bred -= 1.0;
                prod2 *= bred / (ared + bred);
            }
            return fm::log(prod1) + fm::log(prod2)
                    + (logGamma(ared) + (logGamma(bred) - logGammaSum(ared, bred)));
        }

        return fm::log(prod1) + logGamma(ared) + logGammaMinusLogGammaSum(ared, b);
    }

    if(a >= 1.0){
        if(b > 2.0){
            if(b < 10.0){
                double prod = 1.0;
                double bred = b;
                while(bred > 2.0){
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                return fm::log(prod) + (logGamma(a) + (logGamma(bred) - logGammaSum(a, bred)));
            }
            return logGamma(a) + logGammaMinusLogGammaSum(a, b);
        }
        return logGamma(a) + logGamma(b) - logGammaSum(a, b);
    }

    // a < 1 的分支需要 Gamma.gamma()。本研究的 logBeta 參數皆由整數計數而來（>= 1），
    // 故不應觸發；觸發即代表前提錯誤。
    throw std::runtime_error("logBeta a<1 branch reached (unexpected for this study)");
}

double regularizedBeta(double x, double a, double b)
{
    if(std::isnan(x) || std::isnan(a) || std::isnan(b) || x < 0 || x > 1 || a <= 0 || b <= 0){
        return std::numeric_limits<double>::quiet_NaN();
    }

    if(x > (a + 1) / (2 + b + a) && 1 - x <= (b + 1) / (2 + b + a)){
        return 1 - regularizedBeta(1 - x, b, a);
    }

    return fm::exp((a * fm::log(x)) + (b * fm::log1p(-x)) - fm::log(a) - logBeta(a, b))
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
