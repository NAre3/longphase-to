// 開發用 harness：由 CP-A7.tsv 直接跑分段，不重掃 BAM。
// **非移植行為**，與 -debug_only_chr、noisefloor_from_cp_a5 同性質。
#include <cstdio>
#include <string>
#include <vector>

#include "../CpDump.h"
#include "../Segmentation.h"

int main(int argc, char **argv)
{
    if(argc < 3){
        std::fprintf(stderr, "usage: segmentation_from_cp_a7 <CP-A7.tsv> <cpdump_dir> [pcf_out]\n");
        return 2;
    }

    std::FILE *in = std::fopen(argv[1], "r");
    if(in == nullptr){ std::perror("fopen"); return 2; }

    char line[512];
    if(std::fgets(line, sizeof(line), in) == nullptr){ return 2; }

    std::vector<amber::AmberBAF> bafs;
    bafs.reserve(750000);

    while(std::fgets(line, sizeof(line), in) != nullptr){
        char chr[64], tumorBaf[64], tumorModBaf[64], normalBaf[64];
        int pos, idx, tumorDepth, normalDepth;
        if(std::sscanf(line, "%63[^\t]\t%d\t%d\t%63[^\t]\t%63[^\t]\t%d\t%63[^\t]\t%d",
                chr, &pos, &idx, tumorBaf, tumorModBaf, &tumorDepth, normalBaf, &normalDepth) != 8){
            continue;
        }
        amber::AmberBAF baf;
        baf.chromosome = chr;
        baf.position = pos;
        baf.tumorBAF = std::strtod(tumorBaf, nullptr);
        baf.tumorDepth = tumorDepth;
        baf.normalBAF = std::strtod(normalBaf, nullptr);
        baf.normalDepth = normalDepth;
        bafs.push_back(std::move(baf));
    }
    std::fclose(in);

    std::fprintf(stderr, "loaded %zu BAFs\n", bafs.size());

    const amber::SegmentationResult result = amber::segmentBafs(bafs);

    std::fprintf(stderr, "totalCount=%d penaltyMode=%s arms=%zu\n",
            result.totalCount, result.penaltyMode.c_str(), result.arms.size());

    amber::CpDump::setDir(argv[2]);
    amber::writeSegmentationCheckpoints(result);

    if(argc >= 4){
        amber::writeSegmentsFile(argv[3], result);
    }

    std::size_t segments = 0;
    for(const amber::ArmSegments &arm : result.arms){ segments += arm.segments.size(); }
    std::fprintf(stderr, "segments=%zu\n", segments);

    return 0;
}
