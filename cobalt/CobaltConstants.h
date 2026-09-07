#ifndef COBALT_CONSTANTS_H
#define COBALT_CONSTANTS_H

// CobaltConstants.java（spec G9）
namespace cobalt {

constexpr int WINDOW_SIZE = 1000;
constexpr int PARTITION_SIZE = 100000000;   // 100_000_000
constexpr double GC_MAPPABLE_THRESHOLD = 0.85;

// Doubles.greaterOrEqual(a, b) = (a - b) > -1.0E-10
// 位元組碼：dsub; ldc2_w -1.0E-10; dcmpl; ifle ...  → **嚴格大於**
constexpr double DOUBLES_EPSILON = 1.0e-10;
inline bool doublesGreaterOrEqual(double a, double b) { return (a - b) > -DOUBLES_EPSILON; }

}

#endif
