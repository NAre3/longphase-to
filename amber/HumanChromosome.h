#ifndef AMBER_HUMANCHROMOSOME_H
#define AMBER_HUMANCHROMOSOME_H

#include <string>

// 對應 hmf-common HumanChromosome。只保留 AMBER tumor-only 路徑用得到的部分。
namespace amber {

// 去掉開頭的 "chr"（不分大小寫），對應 RefGenomeFunctions.stripChrPrefix
std::string stripChrPrefix(const std::string &chromosome);

// HumanChromosome.contains：1-22 或 X 或 Y，其餘（含 MT）皆為 false
bool isHumanChromosome(const std::string &chromosome);

}

#endif
