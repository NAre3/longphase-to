#ifndef AMBER_FASTMATH_H
#define AMBER_FASTMATH_H

// commons-math3 3.6.1 的 FastMath.log / log1p / exp 在 C++ 端的等價實作。
//
// 為什麼不能用 std::log / std::exp：實測 FastMath 與 java.lang.Math 本身就不逐位元相同
// （20 萬組取樣中 log1p 有 13142 組、exp 有 457 組、log 有 21 組不同），
// 而 std:: 版本走的是系統 libm。用 std:: 版本會讓 binomial CDF 出現 1~32 ulp 的偏差
// （61366 組對照中有 215 組），而 AMBER 拿 CDF 去比 0.16 / 0.84 的門檻做布林判定，
// 偏差有機會翻轉個別點的歸屬。
//
// 查表由 amber_v4.3.jar 內的 bytecode 以反射倒出（FastMathTables.inc），
// 不是抄上游原始碼——jar 才是產生參考值的實體。
//
// **編譯必須關閉 FMA 收縮**（-ffp-contract=off）：a*b+c 若被融合成 fma，
// 中間結果不再經過捨入，會與 Java 不同。

namespace amber { namespace fm {

double log(double x);
double log(double x, double *hiPrec);   // hiPrec 為長度 2 的陣列或 nullptr
double log1p(double x);
double exp(double x);

}}

#endif
