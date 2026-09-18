#ifndef PURPLE_PURPLEINPUTADAPTER_H
#define PURPLE_PURPLEINPUTADAPTER_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

#include "../amber/AmberOutput.h"
#include "../amber/Segmentation.h"
#include "../cobalt/CobaltOutput.h"
#include "../cobalt/Segmentation.h"

namespace purple {

// ---- 這個檔案**不屬於 purple_port** ----
//
// purple/Makefile 的 SRC 刻意不含 PurpleInputAdapter.cpp：它相依 amber/ 與 cobalt/
// 的標頭，而那兩者會拉進 htslib，purple_port 既不需要也沒有連結 htslib。
// 只有頂層 Makefile 的整合建置會編這個 translation unit。
// （對照 amber/AmberApplication.o 被排除在主 build 之外的同一種安排。）
//
// ---- 為什麼要在這裡重做四位小數捨入 ----
//
// AMBER 與 COBALT 寫檔時把 double 截到四位小數：
//   amber/AmberOutput.cpp:19   format4      -> snprintf("%.4f")
//   amber/AmberOutput.cpp:30   formatOrZero -> 非有限值寫字面 "0"
//   cobalt/CobaltOutput.cpp:17 decimalFormat4 -> "%.4f" 去尾零，NaN 寫 "NaN"，"-0" 正規化成 "0"
// 而 PURPLE 以 std::stod 讀回（purple/PurpleInput.cpp:92、113-125）。
// 亦即**凍結候選 34b8d97 所驗證的那組輸入，是量化過的值**。
//
// 記憶體交接若直接傳全精度 double，PURPLE 拿到的輸入就與凍結候選不同，
// 結果會漂移，而 RUN-P009 的 174/174 證據也不再適用。
// 因此這裡以「格式化成同一個字串再 stod 回來」重現**完全相同的運算**——
// 不是近似，是同一個 double -> 十進位字串 -> double 的來回，
// 所以兩條路徑位元相同是由構造保證的，不是靠比對觀察到的。
//
// PCF 不需要這道手續：PURPLE 的 PCF 讀取只用 chromosome/start/end 三個整數欄，
// MeanRatio 讀進來就丟掉（PurpleInput.cpp:131-158），故交接無損。

// 對應 AmberBAFFile 的欄位捨入（tumorBAF/normalBAF 走 format4，
// 非有限值走 formatOrZero 的字面 "0"）。
double roundTripAmberBaf(double value);

// 對應 CobaltRatioFile 的 DecimalFormat("#.####")。
double roundTripCobaltRatio(double value);

// 以 AMBER/COBALT 的記憶體結果組出 InputData，等同「寫出五個 stage 檔案再讀回」。
//
// contamination 對應 amber.qc 的 Contamination 欄；tumor-only 流程固定為 0，
// 由呼叫端傳入而不在此硬編碼，以免整合版與 amber_port 的來源分岔。
InputData buildInputsFromMemory(
        const std::string &sampleId,
        const std::vector<amber::AmberBAF> &bafs,
        const amber::SegmentationResult &amberSegmentation,
        const std::vector<cobalt::CobaltRatio> &ratios,
        const cobalt::SegmentationResult &cobaltSegmentation,
        double contamination);

}

#endif
