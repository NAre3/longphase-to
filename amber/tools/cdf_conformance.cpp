// 數值對照：讀入 Java 端產生的 (n, p, k, cdfBits)，以 C++ 實作重算並比對。
// p 以十六進位浮點字面值傳遞，避免十進位轉換造成的誤差混淆比對結果。
//
// 【2026-09-06】原本逐位元相同即為驗收門檻，離開碼在有任何一筆不同時為 1。
// 改用標準庫後逐位元不同是**預期狀態**（15021/61366），該離開碼會把正常建置報成失敗。
// 現行契約：以相對誤差判定，門檻預設為研究規格 §5 的 1e-9，可由第二個引數覆蓋。
// 逐位元的不一致數只作為診斷輸出。
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

#include "../CommonsMath.h"

int main(int argc, char **argv)
{
    if(argc < 2){
        std::fprintf(stderr, "usage: cdf_conformance <java_cdf.tsv> [relative_tolerance]\n");
        return 2;
    }

    std::FILE *in = std::fopen(argv[1], "r");
    if(in == nullptr){
        std::perror("fopen");
        return 2;
    }

    char line[512];
    if(std::fgets(line, sizeof(line), in) == nullptr){ // header
        return 2;
    }

    const double tolerance = (argc >= 3) ? std::strtod(argv[2], nullptr) : 1e-9;

    long total = 0;
    long mismatches = 0;      // 逐位元不同（診斷用）
    long overTolerance = 0;   // 超出相對容差（驗收用）
    long shown = 0;
    double maxRelative = 0.0;

    while(std::fgets(line, sizeof(line), in) != nullptr){
        int n = 0, k = 0;
        char pHex[128];
        char bits[64];
        if(std::sscanf(line, "%d\t%127[^\t]\t%d\t%63[^\t]", &n, pHex, &k, bits) != 4){
            continue;
        }

        const double p = std::strtod(pHex, nullptr);
        const double cdf = amber::cm3::binomialCdf(n, p, k);

        std::uint64_t cppBits;
        std::memcpy(&cppBits, &cdf, sizeof(cppBits));

        const std::uint64_t javaBits = std::strtoull(bits, nullptr, 16);

        ++total;
        if(cppBits != javaBits){
            ++mismatches;

            double javaValue;
            std::memcpy(&javaValue, &javaBits, sizeof(javaValue));
            const double denominator = std::fabs(javaValue) > std::fabs(cdf)
                    ? std::fabs(javaValue) : std::fabs(cdf);
            const double relative = denominator > 0.0
                    ? std::fabs(cdf - javaValue) / denominator : 0.0;
            if(relative > maxRelative){
                maxRelative = relative;
            }
            if(relative > tolerance){
                ++overTolerance;
            }

            if(relative > tolerance && shown < 10){
                ++shown;
                double javaVal;
                std::memcpy(&javaVal, &javaBits, sizeof(javaVal));
                std::fprintf(stderr, "  n=%d p=%s k=%d java=%.17g cpp=%.17g ulpDiff=%lld\n",
                        n, pHex, k, javaVal, cdf,
                        static_cast<long long>(cppBits) - static_cast<long long>(javaBits));
            }
        }
    }

    std::fclose(in);
    std::printf("compared %ld values, bitwise mismatches %ld, max relative error %.3g\n",
            total, mismatches, maxRelative);
    std::printf("over tolerance %.3g: %ld\n", tolerance, overTolerance);
    return overTolerance == 0 ? 0 : 1;
}
