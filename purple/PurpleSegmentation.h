#ifndef PURPLE_SEGMENTATION_H
#define PURPLE_SEGMENTATION_H

#include <string>
#include <unordered_map>
#include <vector>

#include "PurpleTypes.h"

namespace purple {

// 染色體名 -> 長度。只留 lp::isHumanChromosome 為真者。
using ChromosomeLengths = std::unordered_map<std::string, int>;

// purple_port 的來源：-ref_genome 指向的 fasta 的 .fai 索引（**不讀序列本身**，
// 只取前兩欄）。整合版不走這條：longphase-to 已經開著 BAM，染色體長度直接
// 取自 header @SQ（COBALT 的 loadChromosomes 早就這麼做），不必再要一次
// -ref_genome。
ChromosomeLengths readChromosomeLengths(const std::string &referenceFasta);

std::vector<SupportSegment> createSupportSegments(const InputData &data, const ChromosomeLengths &lengths);
void dumpSupportSegments(const std::vector<SupportSegment> &segments);

}

#endif
