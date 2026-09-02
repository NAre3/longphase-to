#ifndef AMBER_AMBEROUTPUT_H
#define AMBER_AMBEROUTPUT_H

#include <string>
#include <vector>

namespace amber {

// 對應 hmf-common 的 AmberBAF（AmberBAF.java）。
struct AmberBAF
{
    std::string chromosome;
    int position = 0;
    double tumorBAF = 0.0;
    int tumorDepth = 0;
    double normalBAF = 0.0;
    int normalDepth = 0;

    // AmberBAF.java:tumorModifiedBAF / normalModifiedBAF
    double tumorModifiedBAF() const;
    double normalModifiedBAF() const;
};

// 對應 AmberBAFFile.write（4 位小數、非有限值寫 "0"）。以 gzip 寫出。
void writeAmberBafFile(const std::string &path, const std::vector<AmberBAF> &bafs);

// 對應 AmberQCFile.write + AmberQC.status()。
void writeAmberQcFile(const std::string &path, double contamination, double consanguinityProportion);

// 對應 PositionEvidenceFile.write（`-write_tumor_data` 時寫出 <sample>.amber.tumor.raw.tsv.gz）。
// 內容即 CP-A5 的 rawData——同一份已逐位元組驗證過的資料，只是換上 AMBER 的欄名。
// **純輸出，不影響任何計算**：產生 Java 參考端時本就帶著 -write_tumor_data，而 C++ 端在
// 未實作此輸出的情況下，30 組樣本的十一個 checkpoint 與三個 stage 輸出仍逐位元組相同。
struct RawTumorRow
{
    const std::string *chromosome;
    int position;
    char ref;
    char alt;
    int readDepth;
    int indelCount;
    int refSupport;
    int altSupport;
    int baseQualFiltered;
    int mapQualFiltered;
    int seqTechFiltered;
};
void writeTumorRawFile(const std::string &path, const std::vector<RawTumorRow> &rows);

// 對應 VersionInfo.write 產生的 amber.version。
// **刻意不寫成 version=4.3**：那會讓這個檔案看起來像 hmftools 的 AMBER 產出。
// 這裡記錄的是「本移植所重現的 AMBER 版本」，並標明實作者。
void writeVersionFile(const std::string &path);

}

#endif
