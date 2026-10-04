#ifndef AMBER_AMBERSITESFILE_H
#define AMBER_AMBERSITESFILE_H

#include <string>
#include <vector>

#include "AmberSite.h"

namespace amber {

// 對應 AmberSitesFile.loadFile。回傳順序即檔案順序（Java 端為 ArrayListMultimap
// 依 chromosome 分組後的插入順序；CP-A1 dump 前兩端都會再排序，故此處只需保持檔案順序穩定）。
std::vector<AmberSite> loadAmberSites(const std::string &filename);

}

#endif
