#ifndef PURPLE_INPUT_H
#define PURPLE_INPUT_H

#include <string>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

std::vector<AmberBaf> readAmberBafs(const std::string &path);
std::vector<CobaltRatio> readCobaltRatios(const std::string &path, Gender gender);
std::vector<PcfPosition> readPcfPositions(const std::string &path, PcfSource source);
double readAmberContamination(const std::string &path);
Gender determineAmberGender38(const std::vector<AmberBaf> &bafs);
int averageTumorDepth(const std::vector<AmberBaf> &bafs);
InputData loadTumorOnlyInputs(const std::string &sampleId, const std::string &amberDir, const std::string &cobaltDir);
void dumpInputCheckpoint(const InputData &data);

}

#endif
