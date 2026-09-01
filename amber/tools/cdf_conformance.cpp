// EXP-006 數值一致性對照：讀入 Java 端產生的 (n, p, k, cdfBits)，以 C++ 實作重算並逐位元比對。
// p 以十六進位浮點字面值傳遞，避免十進位轉換造成的誤差混淆比對結果。
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "../CommonsMath.h"

int main(int argc, char **argv)
{
    if(argc < 2){
        std::fprintf(stderr, "usage: cdf_conformance <java_cdf.tsv>\n");
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

    long total = 0;
    long mismatches = 0;
    long shown = 0;

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
            if(shown < 10){
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
    std::printf("compared %ld values, bitwise mismatches %ld\n", total, mismatches);
    return mismatches == 0 ? 0 : 1;
}
