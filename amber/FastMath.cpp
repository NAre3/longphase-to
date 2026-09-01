#include "FastMath.h"

#include <cmath>
#include <cstdint>
#include <cstring>

namespace amber { namespace fm {

namespace {

#include "FastMathTables.inc"

constexpr double HEX_40000000 = 1073741824.0; // FastMath.java:307，long 0x40000000 在算式中提升為 double
constexpr double TWO_POWER_52 = 4503599627370496.0;
constexpr double LN_2_A = 0.693147063255310059;
constexpr double LN_2_B = 1.17304635250823482e-7;
constexpr double F_1_3 = 1.0 / 3.0;
constexpr double F_1_2 = 1.0 / 2.0;
constexpr int EXP_INT_TABLE_MAX_INDEX = 750;

constexpr double LN_QUICK_COEF[9][2] = {
    {1.0, 5.669184079525E-24},
    {-0.25, -0.25},
    {0.3333333134651184, 1.986821492305628E-8},
    {-0.25, -6.663542893624021E-14},
    {0.19999998807907104, 1.1921056801463227E-8},
    {-0.1666666567325592, -7.800414592973399E-9},
    {0.1428571343421936, 5.650007086920087E-9},
    {-0.12502530217170715, -7.44321345601866E-11},
    {0.11113807559013367, 9.219544613762692E-9},
};

constexpr double LN_HI_PREC_COEF[6][2] = {
    {1.0, -6.032174644509064E-23},
    {-0.25, -0.25},
    {0.3333333134651184, 1.9868161777724352E-8},
    {-0.2499999701976776, -2.957007209750105E-8},
    {0.19999954104423523, 1.5830993332061267E-10},
    {-0.16624879837036133, -2.6033824355191673E-8},
};

inline std::int64_t rawBits(double x)
{
    std::int64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    return bits;
}

inline double fromBits(std::int64_t bits)
{
    double x;
    std::memcpy(&x, &bits, sizeof(x));
    return x;
}

// FastMath.java: private static double exp(double x, double extra, double[] hiPrec)
double expImpl(double x, double extra, double *hiPrec)
{
    double intPartA;
    double intPartB;
    int intVal = static_cast<int>(x);

    if(x < 0.0){
        if(x < -746.0){
            if(hiPrec != nullptr){ hiPrec[0] = 0.0; hiPrec[1] = 0.0; }
            return 0.0;
        }
        if(intVal < -709){
            const double result = expImpl(x + 40.19140625, extra, hiPrec) / 285040095144011776.0;
            if(hiPrec != nullptr){ hiPrec[0] /= 285040095144011776.0; hiPrec[1] /= 285040095144011776.0; }
            return result;
        }
        if(intVal == -709){
            const double result = expImpl(x + 1.494140625, extra, hiPrec) / 4.455505956692756620;
            if(hiPrec != nullptr){ hiPrec[0] /= 4.455505956692756620; hiPrec[1] /= 4.455505956692756620; }
            return result;
        }
        intVal--;
    }else{
        if(intVal > 709){
            if(hiPrec != nullptr){ hiPrec[0] = INF_D; hiPrec[1] = 0.0; }
            return INF_D;
        }
    }

    intPartA = EXP_INT_TABLE_A[EXP_INT_TABLE_MAX_INDEX + intVal];
    intPartB = EXP_INT_TABLE_B[EXP_INT_TABLE_MAX_INDEX + intVal];

    const int intFrac = static_cast<int>((x - intVal) * 1024.0);
    const double fracPartA = EXP_FRAC_TABLE_A[intFrac];
    const double fracPartB = EXP_FRAC_TABLE_B[intFrac];

    const double epsilon = x - (intVal + intFrac / 1024.0);

    double z = 0.04168701738764507;
    z = z * epsilon + 0.1666666505023083;
    z = z * epsilon + 0.5000000000042687;
    z = z * epsilon + 1.0;
    z = z * epsilon + -3.940510424527919E-20;

    const double tempA = intPartA * fracPartA;
    const double tempB = intPartA * fracPartB + intPartB * fracPartA + intPartB * fracPartB;

    const double tempC = tempB + tempA;

    if(tempC == INF_D){
        return INF_D;
    }

    double result;
    if(extra != 0.0){
        result = tempC * extra * z + tempC * extra + tempC * z + tempB + tempA;
    }else{
        result = tempC * z + tempB + tempA;
    }

    if(hiPrec != nullptr){
        hiPrec[0] = tempA;
        hiPrec[1] = tempC * extra * z + tempC * extra + tempC * z + tempB;
    }

    return result;
}

}

// FastMath.java: private static double log(final double x, final double[] hiPrec)
double log(double x, double *hiPrec)
{
    if(x == 0){
        return -INF_D;
    }

    std::int64_t bits = rawBits(x);

    if(((bits & static_cast<std::int64_t>(0x8000000000000000ULL)) != 0 || x != x) && x != 0.0){
        if(hiPrec != nullptr){ hiPrec[0] = NAN_D; }
        return NAN_D;
    }

    if(x == INF_D){
        if(hiPrec != nullptr){ hiPrec[0] = INF_D; }
        return INF_D;
    }

    int exp = static_cast<int>(bits >> 52) - 1023;

    if((bits & 0x7ff0000000000000LL) == 0){
        if(x == 0){
            if(hiPrec != nullptr){ hiPrec[0] = -INF_D; }
            return -INF_D;
        }
        bits <<= 1;
        while((bits & 0x0010000000000000LL) == 0){
            --exp;
            bits <<= 1;
        }
    }

    if((exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hiPrec == nullptr){
        double xa = x - 1.0;
        double xb = xa - x + 1.0;
        double tmp = xa * HEX_40000000;
        double aa = xa + tmp - tmp;
        double ab = xa - aa;
        xa = aa;
        xb = ab;

        double ya = LN_QUICK_COEF[8][0];
        double yb = LN_QUICK_COEF[8][1];
        for(int i = 7; i >= 0; --i){
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;

            aa = ya + LN_QUICK_COEF[i][0];
            ab = yb + LN_QUICK_COEF[i][1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }

        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;
        tmp = aa * HEX_40000000;
        ya = aa + tmp - tmp;
        yb = aa - ya + ab;

        return ya + yb;
    }

    const double *lnm = LN_MANT[static_cast<int>((bits & 0x000ffc0000000000LL) >> 42)];

    const double epsilon = static_cast<double>(bits & 0x3ffffffffffLL)
            / (TWO_POWER_52 + static_cast<double>(bits & 0x000ffc0000000000LL));

    double lnza = 0.0;
    double lnzb = 0.0;

    if(hiPrec != nullptr){
        double tmp = epsilon * HEX_40000000;
        double aa = epsilon + tmp - tmp;
        double ab = epsilon - aa;
        double xa = aa;
        double xb = ab;

        const double numer = static_cast<double>(bits & 0x3ffffffffffLL);
        const double denom = TWO_POWER_52 + static_cast<double>(bits & 0x000ffc0000000000LL);
        aa = numer - xa * denom - xb * denom;
        xb += aa / denom;

        double ya = LN_HI_PREC_COEF[5][0];
        double yb = LN_HI_PREC_COEF[5][1];
        for(int i = 4; i >= 0; --i){
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;

            aa = ya + LN_HI_PREC_COEF[i][0];
            ab = yb + LN_HI_PREC_COEF[i][1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }

        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;

        lnza = aa + ab;
        lnzb = -(lnza - aa - ab);
    }else{
        lnza = -0.16624882440418567;
        lnza = lnza * epsilon + 0.19999954120254515;
        lnza = lnza * epsilon + -0.2499999997677497;
        lnza = lnza * epsilon + 0.3333333333332802;
        lnza = lnza * epsilon + -0.5;
        lnza = lnza * epsilon + 1.0;
        lnza *= epsilon;
    }

    double a = LN_2_A * exp;
    double b = 0.0;

    double c = a + lnm[0];
    double d = -(c - a - lnm[0]);
    a = c;
    b += d;

    c = a + lnza;
    d = -(c - a - lnza);
    a = c;
    b += d;

    c = a + LN_2_B * exp;
    d = -(c - a - LN_2_B * exp);
    a = c;
    b += d;

    c = a + lnm[1];
    d = -(c - a - lnm[1]);
    a = c;
    b += d;

    c = a + lnzb;
    d = -(c - a - lnzb);
    a = c;
    b += d;

    if(hiPrec != nullptr){
        hiPrec[0] = a;
        hiPrec[1] = b;
    }

    return a + b;
}

double log(double x)
{
    return log(x, nullptr);
}

// FastMath.java: public static double log1p(final double x)
double log1p(double x)
{
    if(x == -1){
        return -INF_D;
    }
    if(x == INF_D){
        return INF_D;
    }

    if(x > 1e-6 || x < -1e-6){
        const double xpa = 1 + x;
        const double xpb = -(xpa - 1 - x);
        double hiPrec[2];
        const double lores = log(xpa, hiPrec);
        if(std::isinf(lores)){
            return lores;
        }
        const double fx1 = xpb / xpa;
        const double epsilon = 0.5 * fx1 + 1;
        return epsilon * fx1 + hiPrec[1] + hiPrec[0];
    }

    const double y = (x * F_1_3 - F_1_2) * x + 1;
    return y * x;
}

double exp(double x)
{
    return expImpl(x, 0.0, nullptr);
}

}}
