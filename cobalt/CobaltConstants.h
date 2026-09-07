#ifndef COBALT_CONSTANTS_H
#define COBALT_CONSTANTS_H

// CobaltConstants.java（spec G9）
namespace cobalt {

constexpr int WINDOW_SIZE = 1000;
constexpr int PARTITION_SIZE = 100000000;   // 100_000_000
constexpr int MAX_SPARSE_CONSOLIDATE_DISTANCE = 3000000;   // 3_000_000
constexpr double GC_MAPPABLE_THRESHOLD = 0.85;

// GC_BUCKET_MIN/MAX = calcGcBucket(DEFAULT_GC_RATIO_MIN/MAX) = round(0.24*100) / round(0.68*100)
// RUN-C002B 的 CP-C8-summary 實測為 24 / 68
constexpr double DEFAULT_GC_RATIO_MIN = 0.24;
constexpr double DEFAULT_GC_RATIO_MAX = 0.68;

constexpr int GC_BUCKET_MIN = 24;
constexpr int GC_BUCKET_MAX = 68;

// Doubles.greaterOrEqual(a, b) = (a - b) > -1.0E-10
// 位元組碼：dsub; ldc2_w -1.0E-10; dcmpl; ifle ...  → **嚴格大於**
constexpr double DOUBLES_EPSILON = 1.0e-10;
inline bool doublesGreaterOrEqual(double a, double b) { return (a - b) > -DOUBLES_EPSILON; }

}

#endif
