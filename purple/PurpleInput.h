#ifndef PURPLE_INPUT_H
#define PURPLE_INPUT_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

std::vector<AmberBaf> readAmberBafs(const std::string &path);
// 由已量化的欄位值組出 CobaltRatio：套用 referenceReadDepth == -1 的特例與
// genderAdjusted。讀檔路徑與記憶體交接路徑共用，這是兩者等價的依據。
CobaltRatio makeCobaltRatio(const std::string &chromosome, int position,
        double referenceReadDepth, double referenceGcRatio, double referenceGcDiploidRatio,
        double referenceGcContent, double tumorReadDepth, double tumorGcRatio,
        double tumorGcContent, Gender gender);

std::vector<CobaltRatio> readCobaltRatios(const std::string &path, Gender gender);
// PCF 檔一列代表的區間。start/end 已是換算後的最終值（新格式 end = rawEnd + 1），
// 兩種檔案格式的差異在 readPcfPositions 內就收斂掉。
struct PcfInterval {
    std::string chromosome;
    int start = 0;
    int end = 0;
};

// 由區間清單組出 PcfPosition（配對、排序、相鄰夾擠）。
// 讀檔路徑與記憶體交接路徑共用此函式，這是兩條路徑等價的依據。
std::vector<PcfPosition> buildPcfPositions(const std::vector<PcfInterval> &intervals, PcfSource source);

std::vector<PcfPosition> readPcfPositions(const std::string &path, PcfSource source);
double readAmberContamination(const std::string &path);
Gender determineAmberGender38(const std::vector<AmberBaf> &bafs);
int averageTumorDepth(const std::vector<AmberBaf> &bafs);
InputData loadTumorOnlyInputs(const std::string &sampleId, const std::string &amberDir, const std::string &cobaltDir);
void dumpInputCheckpoint(const InputData &data);

}

#endif
