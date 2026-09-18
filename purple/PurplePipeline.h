#ifndef PURPLE_PURPLEPIPELINE_H
#define PURPLE_PURPLEPIPELINE_H

#include <string>

#include "PurpleTypes.h"

namespace purple {

// PURPLE tumor-only 流程的參數。欄位與 purple_port 的 CLI 一一對應
// （PurpleApplication.cpp 的 usage），預設值與該處相同。
struct PipelineConfig
{
    std::string sampleId;
    std::string amberDir;
    std::string cobaltDir;
    std::string refGenome;
    std::string ensemblDataDir;
    std::string outputDir;
    int threads = 1;
};

// ---- 為什麼這裡沒有 prescan / postscan ----
//
// AmberPipeline.h 與 CobaltPipeline.h 都拆成 prescan/(掃描)/postscan，拆點由
// 「BAM 掃描發生在哪裡」決定（EXP-I02、EXP-I03）。PURPLE **不讀 BAM**，
// 它消費的是 AMBER/COBALT 的 stage 輸出，流程中沒有可共用的掃描，
// 因此沿用那個兩段式形狀只會是東施效顰。
//
// 這裡改拆「輸入從哪裡來」：
//   loadInputs      從 AMBER/COBALT 的五個 stage 檔案載入 —— purple_port 走這條
//   runFromInputs   從已備妥的 InputData 跑到六個核心輸出 —— 兩條路徑共用
//
// 整合版會改以 AMBER/COBALT 的記憶體結果組出 InputData 後直接呼叫
// runFromInputs，省掉檔案來回。**組 InputData 時必須套用與寫檔端相同的
// 四位小數捨入**（CobaltOutput.cpp 的 decimalFormat4、AmberOutput.cpp 的
// "%.4f"），否則 PURPLE 拿到的輸入精度與凍結候選不同，結果會漂。
// 見 PurpleInputAdapter 的說明。
//
// purple_port 的 main 是這兩段的薄包裝，**與整合版呼叫同一組函式**。
// 這是「凍結候選（34b8d97）的行為未因整合而改變」這句話的依據。

InputData loadInputs(const PipelineConfig &cfg);

void runFromInputs(const PipelineConfig &cfg, const InputData &inputs);

// loadInputs + runFromInputs 的便利包裝。
void run(const PipelineConfig &cfg);

}

#endif
