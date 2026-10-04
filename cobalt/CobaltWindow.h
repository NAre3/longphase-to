#ifndef COBALT_COBALTWINDOW_H
#define COBALT_COBALTWINDOW_H

#include <string>
#include <vector>

#include "GcBuckets.h"
#include "GcProfile.h"
#include "ReadDepth.h"
#include "WindowStatuses.h"

namespace cobalt {

// CobaltWindow.java。gcBucket 為 -1 表示 Java 端的 null。
struct CobaltWindow
{
    std::string chromosomeShort;      // Chromosome 物件的短名（CP-C7 的 dump 用）
    int position = 0;
    double readDepth = 0.0;
    double gcContent = 0.0;
    int gcBucket = -1;
    bool isInExcludedRegion = false;
    bool isInTargetRegion = false;
    bool gcReplaced = false;          // §5.2 見證欄位 #3（correctedByReferenceValue 是否真的取代）
};

struct WindowBuildStats
{
    long gcReplacedCount = 0;         // 取代分支觸發次數
    long referenceLookupFailures = 0; // catch(Exception) 的 fallback 次數（覆蓋度稽核）
    long bucketedReadings = 0;        // 實際 addReading 的次數
};

// BamCalculation.addReading 的逐 window 迴圈。
// pails 會被就地更新（bucket 統計），依 CP-C1 的染色體順序、window 昇冪處理。
std::vector<CobaltWindow> buildWindows(const std::vector<DepthReading> &depths,
                                       const WindowStatuses &statuses,
                                       const GcProfileData &gcData,
                                       GcPailsList &pails,
                                       WindowBuildStats &stats);

}

#endif
