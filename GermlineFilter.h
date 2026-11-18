#ifndef GERMLINE_FILTER_H
#define GERMLINE_FILTER_H

// PURPLE 之後、寫 VCF 之前的 germline 後置過濾。
//
// 規則照離線版逐條實作（research_log 2026-11-16 explore_private_germline_fp，E16 README 的參數表）：
//   TDF <= kHighPurityThreshold：E09 HG9GNu
//     - 錨點（輸出 VCF 的 LP_PON、GT 0|1／1|0、PS 有值的 SNV）先用 PURPLE 的 p、段 CN 判偽 het（E08），排除；
//     - 重讀 BAM，算每個 PASS SNV 的逐 hap 四格表（E03 count_hap）；
//     - H（同 hap ALT 比例）可判用 H，否則 G（PURPLE 期望 AF）可判用 G（E06 v2）；
//     - 兩者都不可判、且 gnomAD 同 REF／ALT 的 AC>=2 才移除（E09 §6.2）。
//   TDF > kHighPurityThreshold：H、G 關閉，GT 挑到的 SNV 中 gnomAD AC>=2 者移除（E15）。
// 被判 germline 的紀錄 FILTER 改為 LP_germlineHG，GT 欄改用未 phased 的寫法（ParsingBam.cpp buildDataLine）。

#include <string>
#include <vector>

#include "ParsingBam.h"
#include "PhasingProcess.h"
#include "purple/PurpleTypes.h"

namespace germline_filter {

struct Inputs {
    // PURPLE 是否有結果（沒有時 H／G 不作用）
    bool purpleAvailable = false;
    // PURPLE 最佳解的 purity（未捨入；內部依 purity.tsv 的 "%.4f" 取值）
    double purplePurity = 0;
    // 寫進 cnv.somatic.tsv 的同一份 copy numbers
    const std::vector<purple::PurpleCopyNumber> *copyNumbers = nullptr;
    // 與 convertNonGermlineToSomatic 同一個判斷（SomaticRefinementPolicy.h）
    bool highPurity = false;
};

// 跑完整個過濾，把要改標籤的 SNV 填進 relabel。log 寫到 stderr。
void run(const PhasingParameters &params, SnpParser &snpFile, ChrPhasingResult &chrPhasingResult,
         const Inputs &inputs, GermlineRelabelSet &relabel);

}

#endif
