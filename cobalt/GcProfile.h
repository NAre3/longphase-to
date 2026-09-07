#ifndef COBALT_GCPROFILE_H
#define COBALT_GCPROFILE_H

#include <string>
#include <vector>

namespace cobalt {

// 對應 hmf-common GCProfile（由 GCProfileFactory.fromLine 建構）
struct GcProfile
{
    std::string chromosome;      // GC profile 檔的原字串，帶 chr 前綴
    int start = 0;               // = pos + 1
    int end = 0;                 // = pos + WINDOW_SIZE
    double gcContent = 0.0;
    double nonNPercentage = 0.0;
    double mappablePercentage = 0.0;

    bool isMappable() const;     // Doubles.greaterOrEqual(mappablePercentage, 0.85)
};

// 依 HumanChromosome ordinal 排序的染色體群組；群組內維持檔案順序。
struct GcProfileData
{
    std::vector<std::string> chromosomes;              // 短名，依 ordinal 昇冪
    std::vector<std::vector<GcProfile>> byChromosome;  // 與 chromosomes 同索引

    const std::vector<GcProfile> *find(const std::string &shortName) const;
};

// GCProfileFactory.loadGCContent(WINDOW_SIZE, path)
GcProfileData loadGcProfile(const std::string &path);

}

#endif
